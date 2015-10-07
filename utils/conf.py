# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.

"""
conf.py -- config/job file parser
"""

class Conf(object):
    """config file reader"""
    def __init__(self):
        self._conf = dict()
        pass
    def read(self, fname):
        """return the sections and options in fname."""
        self._conf.clear()
        ##
        import ConfigParser
        #
        cfg = ConfigParser.RawConfigParser()
        cfglst = cfg.read(fname)
        #
        if( len(cfglst) > 0 ):
            pass
        else:
            raise IOError, "Cannot load file: %s"%(fname)
        ##
        sections = cfg.sections()
        for i in sections:
            options = cfg.options( i )
            td = dict()
            #
            for j in options:
                td[j] = cfg.get(i,j)
                pass
            #
            self._conf[i] = td
            #
            pass
        self._conf['sections'] = sections[:]
        ##
        return self._conf
    pass

##
#--eof--#
