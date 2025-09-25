import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from NuRadioReco.detector import detector
import datetime
import logging
from NuRadioReco.framework.parameters import channelParameters as chp
from scipy.optimize import curve_fit

class Event:
    def __init__(self, station, run, triggerTime, eventNumber, triggerType, wf, time,SNR):
        # Initialize data attributes
        self.triggerTime = triggerTime
        self.eventNumber = eventNumber
        self.wf = wf
        self.time = time
        self.triggerType = triggerType
        self.station=station
        self.run=run
        self.fft=[]
        self.freqs=None
        self.SNR=SNR
        self.envelope=[]
        self.enhanced=False
        
    def plot_wf(self, channel, save_path=None):
        """Plot individual channel waveform (unchanged from original)"""
        if self.wf is None or self.time is None:
            raise ValueError("No waveform data loaded. Call load_data() first.")
        
        if channel < 0 or channel >= len(self.wf):
            raise ValueError(f"Invalid channel number. Must be between 0 and {len(self.wf)-1}.")
        
        plt.figure(figsize=(10, 5))
        plt.plot(self.time[channel], self.wf[channel])
        plt.xlabel('Time [ns]')
        plt.ylabel('Voltage [V]')
        plt.title(f'Station {self.station}, Run {self.run}, Channel {channel}')
        plt.grid(True)
        
        if save_path is not None:
            filename = os.path.join(save_path, 
                                  f'station_{self.station}_run_{self.run}_event_{self.eventNumber}_channel_{channel}_waveform.png')
            plt.savefig(filename, bbox_inches='tight', dpi=300)
            plt.close()
        else:
            plt.show()


    def set_fft(self, ch ,enhanced=False):
        # Compute the FFT of the waveform
        if enhanced==True:
            waveform_fft = np.fft.fft(self.wf_enhanced[ch])
        
            # Compute the corresponding frequency bins
            dt = self.time_enhanced[ch][1] - self.time_enhanced[ch][0]
            freqs = np.fft.fftfreq(len(self.wf_enhanced[ch]), dt)
        else:
            waveform_fft = np.fft.fft(self.wf[ch])
        
            # Compute the corresponding frequency bins
            dt = self.time[ch][1] - self.time[ch][0]
            freqs = np.fft.fftfreq(len(self.wf[ch]), dt)
    
        # Only keep the positive frequencies and corresponding FFT values
        pos_mask = freqs > 0
        pos_freqs = freqs[pos_mask]
        fft_amplitudes = np.abs(waveform_fft[pos_mask])
    
        # Find the frequency with the highest amplitude
        max_idx = np.argmax(fft_amplitudes)
        peak_freq = pos_freqs[max_idx]
        peak_ampl = fft_amplitudes[max_idx]
    
        # Store results
        if not hasattr(self, 'fft_dict'):
            self.fft_dict = {}
        self.freqs = pos_freqs  # Shared across channels (assuming same length)
        self.fft_dict[ch] = fft_amplitudes
    
        #print(f"[FFT] Channel {ch}: Peak frequency = {peak_freq:.3f} GHz, Amplitude = {peak_ampl:.3f}")

    
        
    def plot_fft(self, ch, enhanced=False):
        # Compute or retrieve FFT
        #if not hasattr(self, 'fft_dict') or ch not in self.fft_dict:
        if enhanced==False:
            self.set_fft(ch)
        else:
            self.set_fft(ch,enhanced=True)
            
        freqs = self.freqs
        amplitudes = self.fft_dict[ch]
    
        max_idx = np.argmax(amplitudes)
        peak_freq = freqs[max_idx]
        peak_ampl = amplitudes[max_idx]
        #print(f"[FFT] Channel {ch}: Peak frequency = {peak_freq:.3f} GHz, Amplitude = {peak_ampl:.3f}")
    
        plt.subplot(2, 1, 2)
        plt.plot(freqs, amplitudes, label='FFT Magnitude')
        plt.axvline(x=peak_freq, color='r', linestyle='--', label=f'Peak: {peak_freq:.3f} GHz')
        plt.title(f"Run {self.run}, Event {self.eventNumber}, Channel {ch} - Frequency Domain (FFT)")
        plt.ylabel("Amplitude")
        plt.xlabel("Frequency (GHz)")
        plt.xlim(0, 1.0)
        plt.legend()
        plt.tight_layout()
        plt.show()



    def plot_spectrogram(self, ch, nperseg=128, noverlap=None, cmap='viridis', vmin=None, vmax=None):
        """
        Plot a spectrogram (time-frequency representation) of the waveform data.
        
        Args:
            ch (int): Channel number to plot
            nperseg (int): Length of each segment for FFT (default: 128)
            noverlap (int): Number of points to overlap between segments (default: nperseg//2)
            cmap (str): Matplotlib colormap to use (default: 'viridis')
            vmin (float): Minimum value for color scale (in dB)
            vmax (float): Maximum value for color scale (in dB)
        """
        if noverlap is None:
            noverlap = nperseg // 2
    
        # Compute spectrogram
        freqs, times, Sxx = signal.spectrogram(
            self.wf[ch],
            fs=1/(self.time[ch][1] - self.time[ch][0]),  # Sampling frequency in GHz^-1
            nperseg=nperseg,
            noverlap=noverlap,
            scaling='spectrum'
        )
    
        # Convert to dB
        Sxx_dB = 10 * np.log10(Sxx)
    
        # Limit to frequencies up to 1 GHz
        freq_mask = freqs <= 1.0  # GHz
        freqs = freqs[freq_mask]
        Sxx_dB = Sxx_dB[freq_mask, :]
    
        # Create plot
        plt.figure(figsize=(12, 6))
        plt.pcolormesh(times, freqs, Sxx_dB, shading='auto', cmap=cmap, vmin=vmin, vmax=vmax)
        plt.colorbar(label='Power Spectral Density (dB)')
        plt.title(f'Spectrogram - Station {self.station}, Run {self.run}, Event {self.eventNumber}, Channel {ch}')
        plt.ylabel('Frequency [GHz]')
        plt.xlabel('Time [ns]')
        plt.tight_layout()
        plt.show()


    def get_hilbert_envelope(self, plot_all=False, plot_ch=None):

        """
        Compute the Hilbert envelope for all channels' waveforms.
        
        Args:
            plot_all (bool): Whether to plot all envelopes (default: False)
            plot_ch (int or list): Specific channel(s) to plot (overrides plot_all if specified)
            
        Returns:
            dict: {channel_number: (envelope, time)} pairs for all channels
        """
        envelopes = {}
        
        for ch in range(len(self.wf)):
            # Compute analytic signal and envelope
            analytic_signal = signal.hilbert(self.wf[ch])
            envelope = np.abs(analytic_signal)
            envelopes[ch] = envelope
            
            # Handle plotting
            should_plot = plot_all or (isinstance(plot_ch, int) and (ch == plot_ch) or \
                         (isinstance(plot_ch, list) and (ch in plot_ch)))
            
            if should_plot:
                plt.figure(figsize=(10, 5))
                plt.plot(self.time[ch], self.wf[ch], label='Original', alpha=0.5)
                plt.plot(self.time[ch], envelope, label='Envelope', linewidth=1.5)
                plt.xlabel('Time [ns]')
                plt.ylabel('Voltage [V]')
                plt.title(f'Hilbert Envelope - Ch{ch} (Station {self.station}, Run {self.run})')
                plt.legend()
                plt.grid(True)
                plt.show()
                
        #self.envelope=envelopes
        self.envelope = [envelopes[ch]-np.mean(envelopes[ch]) for ch in range(len(self.wf))]
        return envelopes
        
    def cross_correlate_channels(self, ch1, ch2, max_lag=None, plot=False, envelopes=False):
        """
        Compute cross-correlation between two channels and return the time lag of maximum correlation.
        
        Args:
            ch1 (int): First channel number
            ch2 (int): Second channel number
            max_lag (float, optional): Maximum time lag to consider (in ns). If None, uses full waveform length.
            plot (bool): Whether to plot the cross-correlation function
            
        Returns:
            tuple: (time_lag, correlation_coefficient) where:
                - time_lag: Time difference (ch2 - ch1) in ns where max correlation occurs
                - correlation_coefficient: Normalized correlation value at this lag (0-1)
        """
        # Get waveforms and time arrays
        if envelopes:
            wf1 = self.envelope[ch1]
            wf2 = self.envelope[ch2]
        else:
            wf1 = self.wf[ch1]
            wf2 = self.wf[ch2]
        time1 = self.time[ch1]
        time2 = self.time[ch2]
        
        # Check sampling rates are consistent
        dt1 = time1[1] - time1[0]
        dt2 = time2[1] - time2[0]
        t_diff=time2[0]-time1[0]
        #print(t_diff)
        
        if not np.isclose(dt1, dt2):
            raise ValueError(f"Channel {ch1} and {ch2} have different sampling rates: {dt1} ns vs {dt2} ns")
        
        # Compute cross-correlation using scipy
        correlation = signal.correlate(wf1, wf2, mode='full', method='auto')
        lags = signal.correlation_lags(len(wf1), len(wf2), mode='full') * dt1 - t_diff
        
        # Find the lag with maximum correlation
        max_idx = np.argmax(np.abs(correlation))
        max_lag_value = lags[max_idx]
        max_corr_value = correlation[max_idx]
        
        # Normalize correlation to [0,1]
        norm_corr = np.abs(correlation) / np.sqrt(np.sum(wf1**2)) * np.sqrt(np.sum(wf2**2))
        max_norm_corr = norm_corr[max_idx]
        
        if plot:
            plt.figure(figsize=(10, 5))
            plt.plot(lags, norm_corr)
            plt.axvline(max_lag_value, color='r', linestyle='--', 
                       label=f'Max at {max_lag_value:.2f} ns')
            plt.xlabel('Time Lag [ns]')
            plt.ylabel('Normalized Cross-Correlation')
            plt.title(f'Cross-Correlation: Ch{ch1} vs Ch{ch2}\n'
                     f'Max correlation at {max_lag_value:.2f} ns (strength: {max_norm_corr:.2f})')
            plt.legend()
            plt.grid(True)
            plt.show()
        
        return max_lag_value, max_norm_corr

    def signal_enhancer(self, f_sgl=403, plot=False):
        """
        Enhance signal content at a specific frequency using cross-correlation 
        with a sine wave of half the length of the original signal.
    
        Args:
            f_sgl (float): Target frequency in MHz to enhance.
            plot (bool): Whether to plot original and enhanced signal.
    
        Returns:
            list: List of correlation outputs (same length as original waveform) for each channel
        """
        if self.enhanced==True:
            print("Already applied signal_enhancer")
            return
            
        V_A = 0.005
        f_Hz = f_sgl * 1e6  # Convert MHz to Hz
        results = []
        t_results = []
        
        for ch in range(len(self.wf)):
            t_full = self.time[ch] * 1e-9  # ns → s
            signal_data = self.wf[ch]
            N = len(signal_data)
    
            # Generate sine reference of half-length
            half_idx = N // 2
            t_half = t_full[:half_idx]
            T_A = V_A * np.sin(2 * np.pi * f_Hz * t_half)
            # T_A_full = V_A * np.sin(2 * np.pi * f_Hz * t_full)
            
            # Perform full cross-correlation
            corr_full = signal.correlate(signal_data, T_A, mode='full')
    
            # Extract the central part with same length as original signal
            start = len(corr_full)//3-1
            end = 2*len(corr_full)//3
            corr_central = corr_full[start:end]
            t_result=self.time[ch][half_idx-1:half_idx-1+end-start]
            # print(len(corr_full))
            # print(len(self.time[ch]))
            
            results.append(corr_central)
            t_results.append(t_result)

            
            if plot:
                plt.figure(figsize=(12, 6))
                plt.plot(self.time[ch], signal_data, label='Original Signal', alpha=0.5)
                plt.plot(t_result, corr_central, label=f'Enhanced at {f_sgl:.2f} MHz', linewidth=1.5)
                plt.xlabel('Time [ns]')
                plt.ylabel('Amplitude')
                plt.title(f'Signal Enhancement via Cross-Correlation (Channel {ch})')
                plt.legend()
                plt.grid(True)
                plt.tight_layout()
                plt.show()
                
        self.wf_enhanced=results
        self.time_enhanced=t_results
        self.enhanced=True
        
        return results, t_results

    
    # def cross_correlate_across_chs(self, ch1, ch2, plot=False, return_envelope=False):
    #     """
    #     Cross‑correlate channel ``ch1`` waveform with the **first half** of channel
    #     ``ch2``, compute its Hilbert envelope, and *normalise* the correlation by
    #     dividing it by that envelope.

    #     The result is a phase‑preserving, amplitude‑unity cross‑correlation which
    #     facilitates peak picking and thresholding.

    #     If the waveforms have not yet been enhanced via :py:meth:`signal_enhancer`,
    #     the method will call it automatically.

    #     Parameters
    #     ----------
    #     ch1 : int
    #         Index of the first channel.
    #     ch2 : int
    #         Index of the second channel (only the first half is used for the
    #         correlation as per project convention).
    #     plot : bool, optional
    #         When *True*, plots the **normalised** cross‑correlation and its Hilbert
    #         envelope on the same axes.  Default is *False*.
    #     return_envelope : bool, optional
    #         When *True*, the Hilbert envelope is also returned.  Default is
    #         *False*.

    #     Returns
    #     -------
    #     lags : ndarray
    #         Lag values in nanoseconds.
    #     norm_corr : ndarray
    #         Cross‑correlation divided by its Hilbert envelope (unit amplitude).
    #     envelope : ndarray, optional
    #         The Hilbert envelope (only when *return_envelope* is *True*).
    #     """
    #     # ------------------------------------------------------------------
    #     # Ensure we are working with enhanced waveforms
    #     # ------------------------------------------------------------------
    #     if not self.enhanced:
    #         print("[INFO] Waveforms not enhanced. Applying signal_enhancer() …")
    #         self.signal_enhancer()

    #     # ------------------------------------------------------------------
    #     # Prepare data for correlation
    #     # ------------------------------------------------------------------
    #     wf1 = self.wf_enhanced[ch1]
    #     wf2 = self.wf_enhanced[ch2][:len(self.wf_enhanced[ch2]) // 2]  # First half of ch2
        
    #     dt = self.time_enhanced[ch1][1] - self.time_enhanced[ch1][0]
    #     t_diff=self.time_enhanced[ch2][0]-self.time_enhanced[ch1][0]

    #     # ------------------------------------------------------------------
    #     # Cross‑correlation (SciPy handles zero‑padding internally)
    #     # ------------------------------------------------------------------
    #     corr = signal.correlate(wf1, wf2, mode="full")
    #     lags = (
    #         signal.correlation_lags(len(wf1), len(wf2), mode="full") * dt - t_diff
    #     )

    #     # ------------------------------------------------------------------
    #     # Normalise by waveform energy
    #     # ------------------------------------------------------------------
    #     norm_corr = corr / (np.linalg.norm(wf1) * np.linalg.norm(wf2))

    #     # ------------------------------------------------------------------
    #     # Compute Hilbert envelope
    #     # ------------------------------------------------------------------
    #     envelope = np.abs(signal.hilbert(norm_corr))

    #     # ------------------------------------------------------------------
    #     # Divide by envelope to obtain unit‑amplitude correlation
    #     # ------------------------------------------------------------------
    #     norm_corr = np.divide(
    #         norm_corr,
    #         envelope,
    #         out=np.zeros_like(norm_corr),
    #         where=envelope > 0,
    #     )

    #     # ------------------------------------------------------------------
    #     # Optional plotting
    #     # ------------------------------------------------------------------
    #     if plot:
    #         plt.figure(figsize=(10, 5))
    #         plt.plot(lags, norm_corr, label="Normalised X‑corr", alpha=0.8)
    #         #plt.plot(lags, envelope, label="Hilbert envelope", linewidth=1.2)
    #         max_idx = np.argmax(norm_corr)
    #         plt.axvline(
    #             lags[max_idx],
    #             color="r",
    #             linestyle="--",
    #             label=f"Peak @ {lags[max_idx]:.2f} ns",
    #         )
    #         plt.xlabel("Lag [ns]")
    #         plt.ylabel("Amplitude (normalised)")
    #         plt.title(f"Normalised X‑corr & Envelope: ch{ch1} vs first half ch{ch2}")
    #         plt.grid(True, alpha=0.3)
    #         plt.legend()
    #         plt.tight_layout()
    #         plt.show()

    #     # ------------------------------------------------------------------
    #     # Return results
    #     # ------------------------------------------------------------------
    #     if return_envelope:
    #         return lags, norm_corr, envelope
    #     else:
    #         return lags, norm_corr




    # def extract_phase_from_corr(self, lags, norm_corr, f_sgl=403):
    #     """
    #     Extract the phase delay (dt_0) by fitting the cross-correlation
    #     to a cosine function c(t) = cos(2π f_sgl (t - dt_0)).
    #     0 < dt_0 < 2.48ns
    #     Args:
    #         lags (array): Time lag array [ns]
    #         norm_corr (array): Cross-correlation array
    #         f_sgl (float): Balloon signal frequency in MHz (default: 403)
    
    #     Returns:
    #         float: Phase delay dt_0 in ns, constrained to 0 <= dt_0 < (1/f_sgl)
    #     """
        
    #     # ---- 1. Cut to middle portion ----
    #     start = len(norm_corr) // 3 
    #     end = 2 * len(norm_corr) // 3
    #     lags_cut = lags[start:end]
    #     corr_cut = norm_corr[start:end]
    #     print(len(corr_cut))
        
    #     # ---- 1.5 Check length ----
    #     if len(corr_cut) != 1600:
    #         raise ValueError(f"[ERROR] Expected 1600 points in corr_cut, but got {len(corr_cut)}.")
    
    #     # ---- 2. Define model function ----
    #     f_GHz = f_sgl / 1000.0  # MHz → GHz (since lags are in ns)
        
    #     def cos_model(t, dt_0):
    #         return np.cos(2 * np.pi * f_GHz * (t - dt_0))
        
    #     # ---- 3. Fit the model to the data ----
    #     popt, _ = curve_fit(cos_model, lags_cut, corr_cut, p0=[0])
    #     dt_0_fit = popt[0]
        
    #     # ---- 4. Normalize dt_0 to 0 <= dt_0 < 1/f_sgl ----
    #     period_ns = 1.0 / f_GHz  # ns
    #     dt_0_fit = dt_0_fit % period_ns
        
    #     return dt_0_fit


    def balloon_phase_extracter(self, ch1, ch2, f_sgl=403, plot=False):
        """
        Cross-correlate channel ch1 with the first half of ch2, normalize by
        the Hilbert envelope, and extract phase delay dt_0 via cosine fitting.
    
        Parameters
        ----------
        ch1 : int
            Index of the first channel.
        ch2 : int
            Index of the second channel (first half is used).
        f_sgl : float, optional
            Balloon signal frequency in MHz (default: 403).
        plot : bool, optional
            Plot normalized cross-correlation, Hilbert envelope, and cosine fit.
        return_envelope : bool, optional
            Return envelope in addition to lags, norm_corr, and dt_0.
    
        Returns
        -------
        lags : ndarray
            Lag values in nanoseconds.
        norm_corr : ndarray
            Normalized cross-correlation (divided by its Hilbert envelope).
        dt_0_fit : float
            Phase delay in ns (0 <= dt_0 < 1/f_sgl (2.48ns)).
        envelope : ndarray, optional
            Hilbert envelope (only returned if return_envelope=True).
        """
        # ------------------------------------------------------------------
        # Ensure waveforms are enhanced
        # ------------------------------------------------------------------
        if not self.enhanced:
            print("[INFO] Waveforms not enhanced. Applying signal_enhancer() …")
            self.signal_enhancer()
    
        # ------------------------------------------------------------------
        # Prepare data and compute cross-correlation
        # ------------------------------------------------------------------
        wf1 = self.wf_enhanced[ch1]
        wf2 = self.wf_enhanced[ch2][: len(self.wf_enhanced[ch2]) // 2]
        dt = self.time_enhanced[ch1][1] - self.time_enhanced[ch1][0]
        t_diff = self.time_enhanced[ch2][0] - self.time_enhanced[ch1][0]
    
        corr = signal.correlate(wf1, wf2, mode="full")
        lags = signal.correlation_lags(len(wf1), len(wf2), mode="full") * dt - t_diff
    
        norm_corr = corr / (np.linalg.norm(wf1) * np.linalg.norm(wf2))
        envelope = np.abs(signal.hilbert(norm_corr))
    
        # Avoid divide-by-zero
        norm_corr = np.divide(norm_corr, envelope, out=np.zeros_like(norm_corr), where=envelope > 0)
    
        # ------------------------------------------------------------------
        # Extract phase dt_0
        # ------------------------------------------------------------------
        start = len(norm_corr) // 3
        end = 2 * len(norm_corr) // 3
        lags_cut = lags[start:end]
        corr_cut = norm_corr[start:end]
    
        #print(f"[INFO] Middle portion length: {len(corr_cut)}")
        if len(corr_cut) != 1600:
            raise ValueError(f"[ERROR] Expected 1600 points in corr_cut, got {len(corr_cut)}.")
    
        f_GHz = f_sgl / 1000.0  # MHz → GHz
    
        def cos_model(t, dt_0):
            return np.cos(2 * np.pi * f_GHz * (t - dt_0))
    
        popt, _ = curve_fit(cos_model, lags_cut, corr_cut, p0=[0])
        dt_0_fit = popt[0]
    
        period_ns = 1.0 / f_GHz
        dt_0_fit = dt_0_fit % period_ns
    
        # ------------------------------------------------------------------
        # Optional plotting
        # ------------------------------------------------------------------
        if plot:
            plt.figure(figsize=(10, 6))
            plt.plot(lags, norm_corr, label="Normalized X-corr", alpha=0.7)
            plt.plot(lags, envelope, label="Hilbert envelope", alpha=0.5)
            plt.plot(lags_cut, cos_model(lags_cut, dt_0_fit), "r--", label=f"Cosine fit (dt₀={dt_0_fit:.3f} ns)")
            max_idx = np.argmax(norm_corr)
            plt.axvline(lags[max_idx], color="g", linestyle="--", label=f"Peak @ {lags[max_idx]:.2f} ns")
            plt.xlabel("Lag [ns]")
            plt.ylabel("Amplitude (normalized)")
            plt.title(f"X-corr, Envelope & Phase Fit: ch{ch1} vs first half ch{ch2}")
            plt.legend()
            plt.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.show()
    
        # ------------------------------------------------------------------
        # Return results
        # ------------------------------------------------------------------
    
        return dt_0_fit, lags, norm_corr
