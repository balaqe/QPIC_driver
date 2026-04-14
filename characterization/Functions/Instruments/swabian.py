#%% Imports 
from Functions.BasicFunctions import *
os.chdir(r'C:\Program Files\Swabian Instruments\Time Tagger\driver\python')
import TimeTagger
os.chdir(r'\\ait-pdfs.win.dtu.dk\Services\ELEC\comm\QuantumLab\chip_exps\Lab_Codes')
import matplotlib.pyplot as plt

class SwabianTimeTagger():
    
    def __init__(self, channels=None, channels_idler =None, channels_signal =None, delays_ps=None):
        
        self.channels = channels 
        self.channels_idler = channels_idler
        self.channels_signal = channels_signal
        self.delays_ps = delays_ps

        self.tagger = None

        #self.connect()
        # self.connect_server()
        # self.set_delay(self.channels, self.delays_ps)
        # self.pairs_list = None
        # self.generate_channel_pairs()


    # Connect to Swabian Time Tagger via USB
    def connect(self):
        try: 
            self.tagger = TimeTagger.createTimeTagger()
            print("Connected to Swabian Time Tagger {}".format(self.tagger.getSerial()))
            # self.tagger.reset()
            return
        except Exception as e:
            print("Failed to connect to Time Tagger: ", e)
            return 

    # Connect to Swabian Time Tagger over network
    def connect_server(self):
        # print("Connecting to Swabian Time Tagger over network...")
        # return
        print("Search for Time Taggers on the network...")
        servers = TimeTagger.scanTimeTaggerServers()
        print("{} servers found.".format(len(servers)))
        print(servers)
        print('Information about Time Tagger server on localhost:')
        try:
            server_info = TimeTagger.getTimeTaggerServerInfo(servers[0])
            print(server_info)
        except RuntimeError:
            raise Exception('No Time Tagger server available on "localhost" and the default port 41101.')
        print('Connecting to the server on ' + servers[0])
        # Create a TimeTaggerNetwork instance and connect to the server
        self.tagger = TimeTagger.createTimeTaggerNetwork(servers[0])

    # Set software delays for specified channels
    def set_delay(self, channel_delay=None, delay_ps=None):
        # Use instance defaults if no arguments are provided
        if channel_delay is None:
            channel_delay = self.channels
        if delay_ps is None:
            delay_ps = self.delays_ps

        # Check if we actually have data to process
        if channel_delay is not None and delay_ps is not None:
            if len(channel_delay) != len(delay_ps):
                raise ValueError("Channel and delay lists must be the same length.")
            for ch, dt in zip(channel_delay, delay_ps):
                self.tagger.setInputDelay(ch, int(dt))
                print(f"Channel {ch} delay: {dt} ps")
        else:
            print("No channels or delays provided to set_delay.")

    # Given input channels Signals and input channels Idlers, generate all the possible combination (signal, idler)
    def generate_channel_pairs(self):
        if self.channels_signal is None or self.channels_idler is None:
            return
        pairs = []
        for order in range(1, len(self.channels_signal)+1):
            for i in range(len(self.channels_signal)-order+1):
                channel_signal = self.channels_signal[i:i+order]
                for j in range(len(self.channels_idler)-order+1):
                    channel_idler = self.channels_idler[j:j+order]
                    pairs.append(tuple(sorted(list(channel_signal) + list(channel_idler))))
        self.pairs_list = pairs
        print("Generated channel pairs: ", self.pairs_list)

    # Record single counts for specified channels 
    def record_counts(self, channels=None, integration_ms=500):
            if channels is None:
                channels = self.channels
            binwidth_ps = int(integration_ms * 1e9) 
            counter = TimeTagger.Counter(self.tagger, channels, binwidth=binwidth_ps, n_values=1)
            try:
                print(f"Recording rates (Integration: {integration_ms}ms)")
                while True:
                    counter.startFor(binwidth_ps, clear=True)
                    counter.waitUntilFinished()
                    counts = counter.getData() 
                    rates_khz = (counts[:, 0] / (integration_ms / 1000.0)) / 1000.0
                    output = [f"Ch{ch}: {r:,.3f} kHz" for ch, r in zip(channels, rates_khz)]
                    print(" | ".join(output), end="\t\r")
            # This block triggers when you click the Stop button
            except KeyboardInterrupt:
                print("\nMeasurement stopped.")
            except Exception as e:
                print(f"\nAn error occurred: {e}")
            finally:
                del counter


    def record_coincidences(self, window_ps=1000, integration_ms=500):
        if self.pairs_list is None:
            raise Exception("Channel pairs not defined. Please define channels_signal and channels_idler first.")
        
        coinc = TimeTagger.Coincidences(self.tagger, self.pairs_list, window_ps)
        coinc_channels = coinc.getChannels()
        
        binwidth_ps = int(integration_ms * 1e9)
        counter = TimeTagger.Counter(self.tagger, coinc_channels, binwidth=binwidth_ps, n_values=1)
        print(f"Recording Coincidences (Window: {window_ps}ps, Integration: {integration_ms}ms)")
        print(f"Monitoring groups: {self.pairs_list}")
        try:
            while True:
                counter.startFor(binwidth_ps, clear=True)
                counter.waitUntilFinished()
                counts = counter.getData()
                rates_hz = counts[:, 0] / (integration_ms / 1000.0)
                
                output = [f"{group}: {r:,.1f} Hz" for group, r in zip(self.pairs_list, rates_hz)]
                print(" | ".join(output), end="\r", flush=True)

        except KeyboardInterrupt:
            print("\n\nCoincidence measurement stopped.")
        except Exception as e:
            print(f"\nAn error occurred: {e}")
        finally:
            del counter
            del coinc
            print("Hardware resources released.")

    # Record and plot the histogram for a channel pair
    def record_histogram0(self, channels, binwidth_ps = 10, num_bins = 1000, duration_s = 4):
        corr= TimeTagger.Correlation(self.tagger, channels[0], channels[1], binwidth_ps, num_bins)
        corr.startFor(duration_s * 1e12)
        corr.waitUntilFinished()
        index = corr.getIndex()
        counts = corr.getData()
        print(counts)

        plt.figure()
        plt.bar(index * binwidth_ps, counts)
        plt.xlabel('Bin')
        plt.ylabel('Counts')
        plt.title('Histogram')
        plt.show()



    def plot_histogram(self, channel_1, channel_2, binwidth_ps=100, range_ps=10000, integration_ms=200, live=True, accumulate=False):
        n_bins = int(range_ps / binwidth_ps)
        corr = TimeTagger.Correlation(self.tagger, channel_1, channel_2, 
                                      binwidth=binwidth_ps, n_bins=n_bins)
        
        # 1. Create the figure and axis
        fig, ax = plt.subplots(figsize=(10, 5))
        
        # --- MAKE TRANSPARENT ---
        fig.patch.set_alpha(0.0) 
        ax.patch.set_alpha(0.0)  

        time_axis = corr.getIndex()
        
        # Initialize an array to store the sum of all counts if accumulate is True
        total_data = np.zeros(len(time_axis))
        
        line, = ax.step(time_axis, total_data, where='mid', color='white', linewidth=1.5)
        
        ax.set_title(f"Live Correlation: Ch{channel_1} vs Ch{channel_2}", color='white' if live else 'black')
        ax.set_xlabel("Delay (ps)", color='gray')
        ax.set_ylabel("Counts", color='gray')
        
        ax.grid(True, linestyle=':', alpha=0.3, color='gray')
        for spine in ax.spines.values():
            spine.set_edgecolor('gray')
        ax.tick_params(colors='gray')

        dh = display(fig, display_id=True)
        plt.close(fig) 

        try:
            while True:
                # Collect fresh data for this cycle
                corr.startFor(int(integration_ms * 1e9), clear=True)
                corr.waitUntilFinished()
                current_data = corr.getData()
                
                if accumulate:
                    # Sum the new data with the previous history
                    total_data += current_data
                    line.set_ydata(total_data)
                    current_max = np.max(total_data)
                else:
                    # Just show the snapshot
                    line.set_ydata(current_data)
                    current_max = np.max(current_data)
                
                # Dynamic Y-scaling
                if current_max > 0:
                    ax.set_ylim(-0.5, current_max * 1.1)
                
                if live:
                    dh.update(fig)
                else:
                    break
                    
        except KeyboardInterrupt:
            print("\nLive Plot Stopped.")
        finally:
            del corr


    def close(self):
        TimeTagger.freeTimeTagger(self.tagger)
        del self.tagger
