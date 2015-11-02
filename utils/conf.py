"""
conf.py -- config/job file parser

Changelog:
    conf.py ; started at 2014-09-24 16:48:23 EDT by Liang Chen

    2015-10-26 12:13:57 EDT bugfix by Liang Chen, allow uppercased key name and rename section name list key with underscores.

    2015-10-26 13:02:52 EDT minor changes.
"""

class Conf(object):
    """
    config file reader
    """
    def __init__(self):
        self._conf = dict()
        pass
    def read(self, fname):
        """
        read the file specified in fname, and return the sections and options in a dictionary.

        IOError will be raised when the operation fails.
        """
        self._conf.clear()
        ##
        import ConfigParser
        #
        cfg = ConfigParser.RawConfigParser()
        # keep key name untouched;
        cfg.optionxform = str
        #
        cfglst = list()
        try:
            cfglst.extend( cfg.read(fname) )
            pass
        except:
            raise
        #
        if not( len(cfglst) > 0 ):
            raise IOError, "Cannot load file: %s"%( fname )
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
        self._conf['_sections_'] = sections[:]
        ##
        return self._conf
    def sections(self):
        """
        return the list of all section names.
        """
        res = list()
        ##
        try:
            res.extend( self._conf['_sections_'] )
            pass
        except KeyError:
            pass
        ##
        return res
    pass
##
#--eof--#
