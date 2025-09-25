import os
import requests
import uproot
import readRNOGDataMattak
import logging
#from NuRadioReco.detector.RNO_G import rnog_detector
from NuRadioReco.detector import detector
import datetime
import numpy as np
from NuRadioReco.framework.parameters import channelParametersRNOG as chp_rnog
from NuRadioReco.framework.parameters import channelParameters as chp  # Keep for standard parameters
import NuRadioReco.modules.RNO_G.channelGlitchDetector
import NuRadioReco.modules.channelSignalReconstructor
import NuRadioReco.modules.channelResampler
import NuRadioReco.modules.channelBandPassFilter
import NuRadioReco.modules.channelCWNotchFilter
#import NuRadioReco.modules.RNO_G.hardwareResponseIncorporator
import NuRadioReco.modules.RNO_G.dataProviderRNOG  # Add this import
from NuRadioReco.utilities import units, logging as nulogging
from event_v2 import Event  # Import the Event class
from rnog_analysis_tools.glitch_unscrambler import glitch_unscrambler

# Mute the specific warning from channelBlockOffsetFitter
logging.getLogger('NuRadioReco.RNO_G.channelBlockOffsetFitter').setLevel(logging.ERROR)
print("filter updated!")

class Sampler:
    def __init__(self, station_num, run_number, file_path=None):
        """
        Initialize the Sampler with station and run information
        
        Args:
            station_num (int): Station number
            run_number (int): Run number
            file_path (str, optional): Path to local file. If provided, will use this instead of downloading.
        """
        self.station_num = station_num
        self.run_number = run_number
        self.file_path = file_path
        self.url = f'https://rno-g.uchicago.edu/data/satellite/station{station_num}/'
        #self.url = f'https://rno-g.uchicago.edu/data/full/root/station{station_num}/'
        
        self.username = 'rno-g'
        self.password = '100PeVvs@SummitStation'
        self.temp_file = 'temp_file.root'
        
    def _get_file_path(self):
        """Determine the file path to use, either local or downloaded"""
        if self.file_path:
            return self.file_path
        else:
            # Download the data file
            response_combined = requests.get(
                f'{self.url}run{self.run_number}/combined.root',
                auth=(self.username, self.password))
            
            # Clean up any existing temp file first
            self._cleanup_temp_file()
                
            with open(self.temp_file, 'wb') as f:
                f.write(response_combined.content)
            
            return self.temp_file

    def load_data_delayed(self):
        """
        Download and process data, returning Event objects
        
        Returns:
            list: List of Event objects containing the data
        """
        # Initialize detector
        #det = detector.Detector(source="rnog_mongo")
        det =detector.Detector(json_filename="/users/PAS2608/youwei/RNO_G/NuRadioMC/NuRadioReco/detector/RNO_G/RNO_season_2024.json")
        #det.update(datetime.datetime(2022, 10, 1))
        det.update(datetime.datetime(2023, 10, 1))
        #det.update(datetime.datetime(2024, 9, 14))
        # Get time delays for all channels

            
        # Get the appropriate file path
        file_to_use = self._get_file_path()
            
        # Define reader configurations
        reader_kwargs = {
            #"select_triggers": "FORCE",  # Pass trigger selection
            # "select_triggers": "LT",  # Pass trigger selection
            #"selectors": lambda eventInfo: (eventInfo.triggerTime - np.floor(eventInfo.triggerTime)) < 0.1  # Pass lambda selector
        }
        
        # Initialize data provider with reader_kwargs
        dataProviderRNOG = NuRadioReco.modules.RNO_G.dataProviderRNOG.dataProviderRNOG()
        dataProviderRNOG.begin(files=[file_to_use], det=det, reader_kwargs=reader_kwargs)
        
        # initialize additional modules
        channelResampler = NuRadioReco.modules.channelResampler.channelResampler()
        channelResampler.begin()
    
        channelBandPassFilter = NuRadioReco.modules.channelBandPassFilter.channelBandPassFilter()
        channelBandPassFilter.begin()
    
        # channelCWNotchFilter = NuRadioReco.modules.channelCWNotchFilter.channelCWNotchFilter()
        # channelCWNotchFilter.begin()
    
        # hardwareResponseIncorporator = NuRadioReco.modules.RNO_G.hardwareResponseIncorporator.hardwareResponseIncorporator()
        # hardwareResponseIncorporator.begin()
    
        channelSignalReconstructor = NuRadioReco.modules.channelSignalReconstructor.channelSignalReconstructor(log_level=logging.WARNING)
        channelSignalReconstructor.begin()
        
        events = []
        i=0

        # Get event information
        reader_dict = dataProviderRNOG.reader.get_events_information(
            keys=["station", "run", "triggerTime", "triggerType", "eventNumber"]
        )

        # Initialize the glitch detector to access its unscramble method
        glitch_detector = NuRadioReco.modules.RNO_G.channelGlitchDetector.channelGlitchDetector()
        glitch_detector.begin()  # Initialize the glitch detector
        glich_counter=0
        SNR_filter_counter=0
        event_counter=0
        Effective_ch=[0,1,2,3,5,6,7,9,10,11,22,23]
        
        try:
            for event_info, event_data in zip(reader_dict.values(), dataProviderRNOG.run()):
                event_counter+=1

                    
                station = event_data.get_station(self.station_num)
                if station is None:
                    continue
        
                # Process each channel - check for glitches and unscramble if needed
                gliched=False
                for ch, channel in enumerate(station.iter_channels()):
                    if (channel.has_parameter(chp_rnog.glitch) and channel.get_parameter(chp_rnog.glitch)):
                    # if (channel.has_parameter(chp_rnog.glitch) and channel.get_parameter(chp_rnog.glitch)) and (ch in Effective_ch):
                        gliched=True
                        # print(event_data.get_id())
                        #print(ch,'gliched')
                        continue
                        # original_trace = channel.get_trace()
                        # sampling_rate = channel.get_sampling_rate()
                        # unscrambled_trace = glitch_unscrambler.unscramble(original_trace)
                        # channel.set_trace(unscrambled_trace, sampling_rate)
                        
                        #channel.set_parameter(chp_rnog.is_unscrambled, True)
                if gliched:
                    glich_counter+=1
                    continue

                i=i+1
                # # Update detector with current event time
                # det.update(station.get_station_time())s
                # Process the event through the processing chain
                channelResampler.run(event_data, station, det, sampling_rate=10 * units.GHz)
                
                # channelBandPassFilter.run(
                #     event_data, station, det,
                #     passband=[0.35 * units.GHz, 0.45 * units.GHz],
                #     filter_type='butter', order=2)
                
                channelBandPassFilter.run(
                    event_data, station, det,
                    passband=[0.38 * units.GHz, 0.42 * units.GHz],
                    filter_type='butter', order=2)
                
                # hardwareResponseIncorporator.run(event_data, station, det, sim_to_data=False, mode='phase_only')
                
                # channelCWNotchFilter.run(event_data, station, det)
                
                channelSignalReconstructor.run(event_data, station, det)

                
                # Calculate SNR array
                SNR_arr = []
                for channel in station.iter_channels():
                    SNR_arr.append(channel[chp.SNR]['peak_2_peak_amplitude'])

                #SNR filter Requirement:
                # if np.min(SNR_arr[7])<4:
                #     SNR_filter_counter+=1
                #     continue
                
                # #SNR filter Requirement:
                # if np.min(SNR_arr[6:8])<4.5 or np.max(SNR_arr)<5:
                #     SNR_filter_counter+=1
                #     continue

                    
                # Get waveforms and apply time delays
                time,waveform = event_data.get_waveforms()
                waveform=[ch_wf-np.mean(ch_wf) for ch_wf in waveform]

                
                # Create Event object
                events.append(Event(
                    station=event_info['station'],
                    run=event_info['run'],
                    triggerTime=event_info['triggerTime'],
                    eventNumber=int(event_info['eventNumber']),
                    triggerType=event_info['triggerType'],
                    wf=waveform if len(waveform) > 1 else None,
                    time=time if len(time) > 0 else None,
                    SNR=SNR_arr
                ))
                
                # if len(events)>5:
                #     break
                
        finally:
            # Ensure cleanup even if something goes wrong
            print(f"{event_counter} total event processed")
            print(f"filter {glich_counter} events due to gliched")
            print(f" {SNR_filter_counter} event filtered out due to low SNR")
            print(f" {len(events)} events after filtering")
            dataProviderRNOG.end()
            # Only clean up if we downloaded a temporary file
            if not self.file_path:
                self._cleanup_temp_file()
            
        return events

    def _cleanup_temp_file(self):
        """Remove the temporary file if it exists"""
        if os.path.exists(self.temp_file):
            try:
                os.remove(self.temp_file)
            except Exception as e:
                print(f"Warning: Could not delete temp file {self.temp_file}: {e}")