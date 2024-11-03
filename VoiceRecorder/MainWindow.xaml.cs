using System;
using System.IO;
using System.Numerics;
using System.Windows;
using LiveCharts;
using LiveCharts.Defaults;
using LiveCharts.Wpf;
using NAudio.Wave;
using System.Windows.Threading;
using System.Collections.Generic;
using System.Threading.Tasks;
using MathNet.Numerics.IntegralTransforms;
using System.Diagnostics;

namespace VoiceRecorder
{
    public partial class MainWindow : Window
    {
        // NAudio Variables
        private WaveInEvent? waveIn = null;
        private WaveFileWriter? waveWriter = null;
        private string outputFilePath = "";
        private WaveOutEvent? waveOut = null;
        private AudioFileReader? audioFileReader = null;

        // LiveCharts Variables
        public ChartValues<ObservablePoint> WaveformValues { get; set; }
        public ChartValues<ObservablePoint> FrequencyValues { get; set; }
        private readonly DispatcherTimer timer;

        // FFT Variables
        private readonly int fftLength = 1024; // Must be a power of 2
        private readonly double[] fftBuffer;
        private List<double> sampleBuffer = new List<double>();

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

            // Initialize and start the dispatcher timer
            timer = new DispatcherTimer
            {
                Interval = TimeSpan.FromMilliseconds(100) // Update every 100 ms
            };
            timer.Tick += Timer_Tick;
            timer.Start();
        }

        private string GetUniqueFilePath()
        {
            string timestamp = DateTime.Now.ToString("yyyyMMdd_HHmmss");
            return $"recordedAudio_{timestamp}.wav";
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
                waveIn.RecordingStopped += OnRecordingStopped; // Ensure this line is before StartRecording

                // Initialize WaveFileWriter with unique file path
                outputFilePath = GetUniqueFilePath();
                waveWriter = new WaveFileWriter(outputFilePath, waveIn.WaveFormat);

                waveIn.StartRecording();

                // Update UI
                RecordButton.IsEnabled = false;
                StopButton.IsEnabled = true;
                PlayButton.IsEnabled = false;
                SaveButton.IsEnabled = false;

                Debug.WriteLine("Recording started.");
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

                Debug.WriteLine("Recording stopped.");
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
                if (string.IsNullOrEmpty(outputFilePath) || !File.Exists(outputFilePath))
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

                Debug.WriteLine("Playback started.");
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
                if (string.IsNullOrEmpty(outputFilePath) || !File.Exists(outputFilePath))
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
                    Debug.WriteLine($"Audio saved to {destinationPath}");
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
            int bytesPerSample = waveIn!.WaveFormat.BitsPerSample / 8;
            int sampleCount = e.BytesRecorded / bytesPerSample;

            for (int index = 0; index < sampleCount; index++)
            {
                // Assuming 16-bit audio
                short sampleShort = BitConverter.ToInt16(e.Buffer, index * bytesPerSample);
                float sample = sampleShort / 32768f; // Normalize to -1.0 to 1.0

                // Add to waveform data
                Dispatcher.BeginInvoke(new Action(() =>
                {
                    double time = WaveformValues.Count * 0.1; // Adjust based on actual sampling rate
                    WaveformValues.Add(new ObservablePoint(time, sample));
                    if (WaveformValues.Count > 1000) // Limit number of points to prevent memory issues
                        WaveformValues.RemoveAt(0);
                }));

                // Accumulate samples for FFT
                sampleBuffer.Add(sample);

                if (sampleBuffer.Count >= fftLength)
                {
                    // Copy to fftBuffer
                    for (int i = 0; i < fftLength; i++)
                    {
                        fftBuffer[i] = sampleBuffer[i];
                    }

                    // Clear the buffer for next set of samples
                    sampleBuffer.Clear();
                }
            }

            // Calculate RMS for Decibel
            double rms = CalculateRMS(fftBuffer);
            double decibel = RMSToDecibel(rms);
            Dispatcher.BeginInvoke(new Action(() =>
            {
                DecibelProgressBar.Value = decibel;
            }));
        }

        private void OnRecordingStopped(object sender, StoppedEventArgs e)
        {
            try
            {
                waveIn?.Dispose();
                waveIn = null;

                waveWriter?.Dispose();
                waveWriter = null;

                if (e.Exception != null)
                {
                    MessageBox.Show($"Error during recording: {e.Exception.Message}", "Recording Error", MessageBoxButton.OK, MessageBoxImage.Error);
                    Debug.WriteLine($"Recording stopped with error: {e.Exception.Message}");
                }
                else
                {
                    // Enable PlayButton and SaveButton only after recording has stopped successfully
                    Dispatcher.BeginInvoke(new Action(() =>
                    {
                        PlayButton.IsEnabled = true;
                        SaveButton.IsEnabled = true;
                    }));
                    Debug.WriteLine("Recording stopped successfully.");
                }
            }
            catch (Exception ex)
            {
                MessageBox.Show($"Exception in OnRecordingStopped: {ex.Message}", "Exception Error", MessageBoxButton.OK, MessageBoxImage.Error);
                Debug.WriteLine($"Exception in OnRecordingStopped: {ex.Message}");
            }
        }

        private void OnPlaybackStopped(object sender, StoppedEventArgs e)
        {
            try
            {
                waveOut?.Dispose();
                waveOut = null;

                audioFileReader?.Dispose();
                audioFileReader = null;

                // Update UI
                Dispatcher.BeginInvoke(new Action(() =>
                {
                    PlayButton.IsEnabled = true;
                    RecordButton.IsEnabled = true;
                    StopButton.IsEnabled = false;
                    SaveButton.IsEnabled = true;
                }));

                if (e.Exception != null)
                {
                    MessageBox.Show($"Error during playback: {e.Exception.Message}", "Playback Error", MessageBoxButton.OK, MessageBoxImage.Error);
                    Debug.WriteLine($"Playback stopped with error: {e.Exception.Message}");
                }
                else
                {
                    Debug.WriteLine("Playback stopped successfully.");
                }
            }
            catch (Exception ex)
            {
                MessageBox.Show($"Exception in OnPlaybackStopped: {ex.Message}", "Exception Error", MessageBoxButton.OK, MessageBoxImage.Error);
                Debug.WriteLine($"Exception in OnPlaybackStopped: {ex.Message}");
            }
        }

        private async void Timer_Tick(object? sender, EventArgs e)
        {
            if (waveIn == null || fftBuffer.Length == 0) return; // Only process if recording and buffer is filled

            Debug.WriteLine("Performing FFT...");

            // Perform FFT on the buffer asynchronously to prevent UI blocking
            await Task.Run(() =>
            {
                // Prepare a complex buffer for FFT
                Complex[] complexBuffer = new Complex[fftLength];
                for (int i = 0; i < fftLength; i++)
                {
                    complexBuffer[i] = new Complex(fftBuffer[i], 0);
                }

                // Apply a windowing function (e.g., Hamming window) to reduce spectral leakage
                for (int i = 0; i < fftLength; i++)
                {
                    double window = 0.54 - 0.46 * Math.Cos(2 * Math.PI * i / (fftLength - 1)); // Hamming window
                    complexBuffer[i] *= window;
                }

                // Perform FFT using MathNet.Numerics
                Fourier.Forward(complexBuffer, FourierOptions.Matlab);

                // Calculate magnitude
                double[] magnitudes = new double[fftLength / 2];
                for (int i = 0; i < magnitudes.Length; i++)
                {
                    magnitudes[i] = complexBuffer[i].Magnitude;
                }

                // Find the index of the maximum magnitude
                int maxIndex = 0;
                double maxMagnitude = 0;
                for (int i = 0; i < magnitudes.Length; i++)
                {
                    if (magnitudes[i] > maxMagnitude)
                    {
                        maxMagnitude = magnitudes[i];
                        maxIndex = i;
                    }
                }

                // Calculate the dominant frequency
                double dominantFrequency = maxIndex * (waveIn.WaveFormat.SampleRate / (double)fftLength);

                // Optional: Apply parabolic interpolation for better accuracy
                if (maxIndex > 0 && maxIndex < magnitudes.Length - 1)
                {
                    double alpha = magnitudes[maxIndex - 1];
                    double beta = magnitudes[maxIndex];
                    double gamma = magnitudes[maxIndex + 1];

                    // Parabolic interpolation formula
                    double p = 0.5 * (alpha - gamma) / (alpha - 2 * beta + gamma);
                    dominantFrequency = (maxIndex + p) * (waveIn.WaveFormat.SampleRate / (double)fftLength);
                }

                // Prepare data for UI thread
                List<ObservablePoint> frequencyData = new List<ObservablePoint>();
                for (int i = 0; i < magnitudes.Length; i++)
                {
                    double frequency = i * (waveIn.WaveFormat.SampleRate / (double)fftLength);
                    frequencyData.Add(new ObservablePoint(frequency, magnitudes[i]));
                }

                // Update FrequencyChart and CurrentFrequencyLabel on UI thread
                Dispatcher.BeginInvoke(new Action(() =>
                {
                    FrequencyValues.Clear();
                    foreach (var point in frequencyData)
                    {
                        FrequencyValues.Add(point);
                    }

                    // Update the Current Frequency Label
                    CurrentFrequencyLabel.Text = $"{dominantFrequency:F1} Hz";
                }));

                Debug.WriteLine($"FFT completed. Dominant Frequency: {dominantFrequency:F1} Hz");
            });
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

        protected override void OnClosed(EventArgs e)
        {
            base.OnClosed(e);

            waveIn?.Dispose();
            waveWriter?.Dispose();
            waveOut?.Dispose();
            audioFileReader?.Dispose();
            timer.Stop();

            Debug.WriteLine("Application closed and resources disposed.");
        }
    }
}
