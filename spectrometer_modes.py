class mode_800(object):

    def __init__(self):

        self.bitcode = 'adc5g_800mhz_qbs_4inbb_2017_May_26_1449.bof.gz'
        self.clk = 1600
        self.bandwidth = 800
        self.ADCstreams = 8
        self.numchannels = 2048
        self.gain = 0xffffff
        self.shift = 0xffff

        self.sync_LCM = 10
        self.pfb_taps = 4
        self.FFTsize = self.numchannels*2
        self.FFTinputs = 8
        self.nbram = 4

class mode_400(object):

    def __init__(self):

        self.bitcode = 'adc5g_400mhz_qbs_4in_2017_May_16_1346.bof.gz'
        self.clk = 800
        self.bandwidth = 400
        self.ADCstreams = 8
        self.numchannels = 4096
        self.gain = 0xffffff
        self.shift = 0xffff
        
        self.sync_LCM = 10
        self.pfb_taps = 4
        self.FFTsize = self.numchannels*2
        self.FFTinputs = 8
        self.nbram = 4

class mode_200(object):

    def __init__(self):

        self.bitcode = 'adc5g_200mhz_qbs_4in_2017_May_22_1352.bof.gz'
        self.clk = 800
        self.bandwidth = 200
        self.ADCstreams = 4
        self.numchannels = 8192
        self.gain = 0xffffff
        self.shift = 0xffff

        self.sync_LCM = 12
        self.pfb_taps = 4
        self.FFTsize = self.numchannels*2
        self.FFTinputs = 4
        self.nbram = 2
