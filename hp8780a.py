#!/usr/bin/env python

from myGpib import Gpib

class HP8780A(Gpib):
    def __init__(self, name='8780a',
                 pad=None, sad=0):
        Gpib.__init__(self, name=name, pad=pad, sad=sad)
    
    def set_freq(self, freq):
        """Sets CW frequency in Hz"""
        if freq < 1e9:
            fstr = "FR%sMZ" % (freq/1.e6)
        else:
            fstr = "FR%sGZ" % (freq/1.e9)
        self.write(fstr)

    def output_off(self):
        self.write('RF0')

    def output_on(self):
        self.write('RF1')

    def set_power_level(self, power):
        self.write("PL%sDBM" % power)
