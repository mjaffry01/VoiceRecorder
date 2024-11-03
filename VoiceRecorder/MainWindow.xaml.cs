using System;
using System.IO;
using System.Numerics;
using System.Windows;
using LiveCharts;
using LiveCharts.Defaults;
using LiveCharts.Wpf;
using NAudio.Wave;
using System.Windows.Threading;

namespace VoiceRecorder
{
    public partial class MainWindow : Window
    {
        // NAudio Variables
        private WaveInEvent? waveIn = null;
        private WaveFileWriter? waveWriter = null;
        private readonly string outputFilePath = "recordedAudio.wav";
        private WaveOutEvent? waveOut = null;
        private AudioFileReader? audioFileReader = null;

        // LiveCharts Variables
        public ChartValues<ObservablePoint> WaveformValues { get; set; }
        public ChartValues<ObservablePoint> FrequencyValues { get; set; }
        private readonly DispatcherTimer timer;

        // FFT Variables
        private readonly int fftLength = 1024; // Must be a power of 2
        private readonly double[] fftBuffer;
        private readonly Complex[] fftComplexBuffer;

        public MainWindow()
        {
            InitializeComponent();
            StopButton.IsEnabled = false;
            PlayButton.IsEnabled = false;
            SaveButton.IsEnabled = false;

            // Initialize LiveCharts data structures
            WaveformValues = new ChartValues<ObservablePoint>();
            FrequencyValues = new ChartValues<ObservablePoint>();

            // Configure Waveform Chart
            WaveformChart.Series = new SeriesCollection
            {
                new LineSeries
                {
                    Title = "Waveform",
                    Values = WaveformValues,
                    PointGeometry = null,
                    StrokeThickness = 1
                }
            };
            WaveformChart.AxisX[0].Title = "Time (ms)";
            WaveformChart.AxisY[0].Title = "Amplitude";

            // Configure Frequency Spectrum Chart
            FrequencyChart.Series = new SeriesCollection
            {
                new ColumnSeries
                {
                    Title = "Frequency",
                    Values = FrequencyValues,
                    DataLabels = false
                }
            };
            FrequencyChart.AxisX[0].Title = "Frequency (Hz)";
            FrequencyChart.AxisY[0].Title = "Amplitude";

            // Initialize FFT buffers
            fftBuffer = new double[fftLength];
            fftComplexBuffer = new Complex[fftLength];

            // Initialize and start the dispatcher timer
            timer = new DispatcherTimer
            {
                Interval = TimeSpan.FromMilliseconds(100) // Update every 100 ms
            };
            timer.Tick += Timer_Tick;
            timer.Start();
        }

        private void RecordButton_Click(object sender, RoutedEventArgs e)
        {
            try
            {
                // Initialize waveIn for recording
                waveIn = new WaveInEvent
                {
                    WaveFormat = new WaveFormat(44100, 1) // 44.1kHz, Mono
                };
                waveIn.DataAvailable += OnDataAvailable;
                waveIn.RecordingStopped += OnRecordingStopped;

                // Initialize WaveFileWriter to write data to a WAV file
                waveWriter = new WaveFileWriter(outputFilePath, waveIn.WaveFormat);

                waveIn.StartRecording();

                // Update UI
                RecordButton.IsEnabled = false;
                StopButton.IsEnabled = true;
                PlayButton.IsEnabled = false;
                SaveButton.IsEnabled = false;
            }
            catch (Exception ex)
            {
                MessageBox.Show($"Error starting recording: {ex.Message}", "Recording Error", MessageBoxButton.OK, MessageBoxImage.Error);
            }
        }

        private void StopButton_Click(object sender, RoutedEventArgs e)
        {
            try
            {
                waveIn?.StopRecording();

                // Update UI
                RecordButton.IsEnabled = true;
                StopButton.IsEnabled = false;
                // PlayButton and SaveButton are enabled in OnRecordingStopped
            }
            catch (Exception ex)
            {
                MessageBox.Show($"Error stopping recording: {ex.Message}", "Stopping Error", MessageBoxButton.OK, MessageBoxImage.Error);
            }
        }

        private void PlayButton_Click(object sender, RoutedEventArgs e)
        {
            try
            {
                if (!File.Exists(outputFilePath))
                {
                    MessageBox.Show("No recording found to play.", "Playback Error", MessageBoxButton.OK, MessageBoxImage.Warning);
                    return;
                }

                if (IsFileLocked(outputFilePath))
                {
                    MessageBox.Show("The audio file is currently in use. Please try again shortly.", "File In Use", MessageBoxButton.OK, MessageBoxImage.Warning);
                    return;
                }

                // Initialize audio playback
                waveOut = new WaveOutEvent();
                audioFileReader = new AudioFileReader(outputFilePath);
                waveOut.Init(audioFileReader);
                waveOut.PlaybackStopped += OnPlaybackStopped;
                waveOut.Play();

                // Update UI
                PlayButton.IsEnabled = false;
                RecordButton.IsEnabled = false;
                StopButton.IsEnabled = false;
                SaveButton.IsEnabled = false;
            }
            catch (Exception ex)
            {
                MessageBox.Show($"Error playing audio: {ex.Message}", "Playback Error", MessageBoxButton.OK, MessageBoxImage.Error);
            }
        }

        private void SaveButton_Click(object sender, RoutedEventArgs e)
        {
            try
            {
                if (!File.Exists(outputFilePath))
                {
                    MessageBox.Show("No recording found to save.", "Save Error", MessageBoxButton.OK, MessageBoxImage.Warning);
                    return;
                }

                if (IsFileLocked(outputFilePath))
                {
                    MessageBox.Show("The audio file is currently in use. Please try again shortly.", "File In Use", MessageBoxButton.OK, MessageBoxImage.Warning);
                    return;
                }

                // Let user choose the save location
                Microsoft.Win32.SaveFileDialog saveFileDialog = new Microsoft.Win32.SaveFileDialog
                {
                    FileName = "recordedAudio",
                    DefaultExt = ".wav",
                    Filter = "WAV files (*.wav)|*.wav"
                };

                bool? result = saveFileDialog.ShowDialog();

                if (result == true)
                {
                    string destinationPath = saveFileDialog.FileName;
                    File.Copy(outputFilePath, destinationPath, overwrite: true);
                    MessageBox.Show($"Audio saved to {destinationPath}", "Save Successful", MessageBoxButton.OK, MessageBoxImage.Information);
                }
            }
            catch (Exception ex)
            {
                MessageBox.Show($"Error saving audio: {ex.Message}", "Save Error", MessageBoxButton.OK, MessageBoxImage.Error);
            }
        }

        private void OnDataAvailable(object sender, WaveInEventArgs e)
        {
            if (waveWriter == null) return;

            waveWriter.Write(e.Buffer, 0, e.BytesRecorded);
            waveWriter.Flush();

            // Convert bytes to float samples
            int bytesPerSample = waveIn.WaveFormat.BitsPerSample / 8;
            int sampleCount = e.BytesRecorded / bytesPerSample;

            for (int index = 0; index < sampleCount; index++)
            {
                // Assuming 16-bit audio
                short sampleShort = BitConverter.ToInt16(e.Buffer, index * bytesPerSample);
                float sample = sampleShort / 32768f; // Normalize to -1.0 to 1.0

                // Add to waveform data
                Dispatcher.Invoke(() =>
                {
                    double time = WaveformValues.Count * 0.1; // 100 ms per sample (adjust as needed)
                    WaveformValues.Add(new ObservablePoint(time, sample));
                    if (WaveformValues.Count > 1000) // Limit number of points to prevent memory issues
                        WaveformValues.RemoveAt(0);
                });

                // Fill FFT buffer
                fftBuffer[index % fftLength] += sample;
            }

            // Calculate RMS for Decibel
            double rms = CalculateRMS(fftBuffer);
            double decibel = RMSToDecibel(rms);
            Dispatcher.Invoke(() =>
            {
                DecibelProgressBar.Value = decibel;
            });
        }

        private void OnRecordingStopped(object sender, StoppedEventArgs e)
        {
            waveIn?.Dispose();
            waveIn = null;

            waveWriter?.Dispose();
            waveWriter = null;

            if (e.Exception != null)
            {
                MessageBox.Show($"Error during recording: {e.Exception.Message}", "Recording Error", MessageBoxButton.OK, MessageBoxImage.Error);
            }
            else
            {
                // Enable PlayButton and SaveButton only after recording has stopped successfully
                Dispatcher.Invoke(() =>
                {
                    PlayButton.IsEnabled = true;
                    SaveButton.IsEnabled = true;
                });
            }
        }

        private void OnPlaybackStopped(object sender, StoppedEventArgs e)
        {
            waveOut?.Dispose();
            waveOut = null;

            audioFileReader?.Dispose();
            audioFileReader = null;

            // Update UI
            Dispatcher.Invoke(() =>
            {
                PlayButton.IsEnabled = true;
                RecordButton.IsEnabled = true;
                StopButton.IsEnabled = false;
                SaveButton.IsEnabled = true;
            });

            if (e.Exception != null)
            {
                MessageBox.Show($"Error during playback: {e.Exception.Message}", "Playback Error", MessageBoxButton.OK, MessageBoxImage.Error);
            }
        }

        private void Timer_Tick(object? sender, EventArgs e)
        {
            if (waveIn == null) return; // Only process if recording

            // Perform FFT on the buffer
            for (int i = 0; i < fftLength; i++)
            {
                fftComplexBuffer[i] = new Complex(fftBuffer[i], 0);
            }

            // Perform FFT
            FastFourierTransform.FFT(true, (int)Math.Log(fftLength, 2.0), fftComplexBuffer);

            // Calculate magnitude
            double[] magnitudes = new double[fftLength / 2];
            for (int i = 0; i < magnitudes.Length; i++)
            {
                magnitudes[i] = fftComplexBuffer[i].Magnitude;
            }

            // Update FrequencyChart
            Dispatcher.Invoke(() =>
            {
                FrequencyValues.Clear();
                for (int i = 0; i < magnitudes.Length; i++)
                {
                    double frequency = i * (waveIn.WaveFormat.SampleRate / (double)fftLength);
                    FrequencyValues.Add(new ObservablePoint(frequency, magnitudes[i]));
                }
            });

            // Reset FFT buffer
            Array.Clear(fftBuffer, 0, fftBuffer.Length);
        }

        // RMS Calculation for Decibel
        private double CalculateRMS(double[] buffer)
        {
            double sum = 0;
            foreach (var sample in buffer)
            {
                sum += sample * sample;
            }
            return Math.Sqrt(sum / buffer.Length);
        }

        // Convert RMS to Decibels
        private double RMSToDecibel(double rms)
        {
            if (rms <= 0.0001)
                return -60; // Minimum dB
            return 20 * Math.Log10(rms);
        }

        // Check if file is locked
        private bool IsFileLocked(string filePath)
        {
            FileStream? stream = null;

            try
            {
                stream = new FileStream(filePath, FileMode.Open, FileAccess.Read, FileShare.None);
            }
            catch (IOException)
            {
                // The file is unavailable because it is:
                // - Still being written to
                // - Or being processed by another thread
                // - Or does not exist (handled elsewhere)
                return true;
            }
            finally
            {
                stream?.Close();
            }

            // File is not locked
            return false;
        }
    }

    // Fast Fourier Transform Implementation
    public static class FastFourierTransform
    {
        public static void FFT(bool forward, int m, Complex[] data)
        {
            int n = 1 << m;
            int i, j, k, n1, n2, a;
            Complex c, t;
            // Bit-reversal
            j = 0;
            n2 = n / 2;
            for (i = 1; i < (n - 1); i++)
            {
                n1 = n2;
                while (j >= n1)
                {
                    j = j - n1;
                    n1 = n1 / 2;
                }
                j = j + n1;

                if (i < j)
                {
                    t = data[i];
                    data[i] = data[j];
                    data[j] = t;
                }
            }

            // FFT
            n1 = 0;
            n2 = 1;

            for (i = 0; i < m; i++)
            {
                n1 = n2;
                n2 = n2 + n2;
                a = 0;

                for (j = 0; j < n1; j++)
                {
                    double cAngle = -2 * Math.PI * a / n2;
                    if (!forward)
                        cAngle = -cAngle;
                    c = new Complex(Math.Cos(cAngle), Math.Sin(cAngle));

                    for (k = j; k < n; k = k + n2)
                    {
                        t = data[k + n1] * c;
                        data[k + n1] = data[k] - t;
                        data[k] = data[k] + t;
                    }
                    a += 1;
                }
            }

            // If inverse FFT, divide by n
            if (!forward)
            {
                for (i = 0; i < n; i++)
                {
                    data[i] /= n;
                }
            }
        }
    }
}
