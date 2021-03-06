# This file was automatically generated by SWIG (http://www.swig.org).
# Version 2.0.11
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.





from sys import version_info
if version_info >= (2,6,0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_rjmcmc', [dirname(__file__)])
        except ImportError:
            import _rjmcmc
            return _rjmcmc
        if fp is not None:
            try:
                _mod = imp.load_module('_rjmcmc', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _rjmcmc = swig_import_helper()
    del swig_import_helper
else:
    import _rjmcmc
del version_info
try:
    _swig_property = property
except NameError:
    pass # Python < 2.2 doesn't have 'property'.
def _swig_setattr_nondynamic(self,class_type,name,value,static=1):
    if (name == "thisown"): return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'SwigPyObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    if (not static):
        self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)

def _swig_setattr(self,class_type,name,value):
    return _swig_setattr_nondynamic(self,class_type,name,value,0)

def _swig_getattr(self,class_type,name):
    if (name == "thisown"): return self.this.own()
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError(name)

def _swig_repr(self):
    try: strthis = "proxy of " + self.this.__repr__()
    except: strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

try:
    _object = object
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0



def rjmcmc_seed(*args):
  return _rjmcmc.rjmcmc_seed(*args)
rjmcmc_seed = _rjmcmc.rjmcmc_seed
class dataset1d(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, dataset1d, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, dataset1d, name)
    __repr__ = _swig_repr
    __swig_setmethods__["d"] = _rjmcmc.dataset1d_d_set
    __swig_getmethods__["d"] = _rjmcmc.dataset1d_d_get
    if _newclass:d = _swig_property(_rjmcmc.dataset1d_d_get, _rjmcmc.dataset1d_d_set)
    def __init__(self, *args): 
        this = _rjmcmc.new_dataset1d(*args)
        try: self.this.append(this)
        except: self.this = this
    def set_xrange(self, *args): return _rjmcmc.dataset1d_set_xrange(self, *args)
    def get_xmin(self, *args): return _rjmcmc.dataset1d_get_xmin(self, *args)
    def get_xmax(self, *args): return _rjmcmc.dataset1d_get_xmax(self, *args)
    def set_yrange(self, *args): return _rjmcmc.dataset1d_set_yrange(self, *args)
    def get_ymin(self, *args): return _rjmcmc.dataset1d_get_ymin(self, *args)
    def get_ymax(self, *args): return _rjmcmc.dataset1d_get_ymax(self, *args)
    def set_lambda_std(self, *args): return _rjmcmc.dataset1d_set_lambda_std(self, *args)
    def get_lambda_std(self, *args): return _rjmcmc.dataset1d_get_lambda_std(self, *args)
    def set_lambda_range(self, *args): return _rjmcmc.dataset1d_set_lambda_range(self, *args)
    def get_lambda_min(self, *args): return _rjmcmc.dataset1d_get_lambda_min(self, *args)
    def get_lambda_max(self, *args): return _rjmcmc.dataset1d_get_lambda_max(self, *args)
    __swig_destroy__ = _rjmcmc.delete_dataset1d
    __del__ = lambda self : None;
dataset1d_swigregister = _rjmcmc.dataset1d_swigregister
dataset1d_swigregister(dataset1d)

class resultset1d(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, resultset1d, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, resultset1d, name)
    __repr__ = _swig_repr
    __swig_setmethods__["r"] = _rjmcmc.resultset1d_r_set
    __swig_getmethods__["r"] = _rjmcmc.resultset1d_r_get
    if _newclass:r = _swig_property(_rjmcmc.resultset1d_r_get, _rjmcmc.resultset1d_r_set)
    __swig_destroy__ = _rjmcmc.delete_resultset1d
    __del__ = lambda self : None;
    def proposed(self, *args): return _rjmcmc.resultset1d_proposed(self, *args)
    def acceptance(self, *args): return _rjmcmc.resultset1d_acceptance(self, *args)
    def partitions(self, *args): return _rjmcmc.resultset1d_partitions(self, *args)
    def order_histogram(self, *args): return _rjmcmc.resultset1d_order_histogram(self, *args)
    def partition_histogram(self, *args): return _rjmcmc.resultset1d_partition_histogram(self, *args)
    def partition_location_histogram(self, *args): return _rjmcmc.resultset1d_partition_location_histogram(self, *args)
    def x(self, *args): return _rjmcmc.resultset1d_x(self, *args)
    def y(self, *args): return _rjmcmc.resultset1d_y(self, *args)
    def mean(self, *args): return _rjmcmc.resultset1d_mean(self, *args)
    def median(self, *args): return _rjmcmc.resultset1d_median(self, *args)
    def mode(self, *args): return _rjmcmc.resultset1d_mode(self, *args)
    def credible_min(self, *args): return _rjmcmc.resultset1d_credible_min(self, *args)
    def credible_max(self, *args): return _rjmcmc.resultset1d_credible_max(self, *args)
    def misfit(self, *args): return _rjmcmc.resultset1d_misfit(self, *args)
    def lambda_history(self, *args): return _rjmcmc.resultset1d_lambda_history(self, *args)
    def histogram(self, *args): return _rjmcmc.resultset1d_histogram(self, *args)
    def __init__(self): 
        """__init__(resultset1d self) -> resultset1d"""
        this = _rjmcmc.new_resultset1d()
        try: self.this.append(this)
        except: self.this = this
resultset1d_swigregister = _rjmcmc.resultset1d_swigregister
resultset1d_swigregister(resultset1d)

class resultset1dfm(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, resultset1dfm, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, resultset1dfm, name)
    __repr__ = _swig_repr
    __swig_setmethods__["r"] = _rjmcmc.resultset1dfm_r_set
    __swig_getmethods__["r"] = _rjmcmc.resultset1dfm_r_get
    if _newclass:r = _swig_property(_rjmcmc.resultset1dfm_r_get, _rjmcmc.resultset1dfm_r_set)
    __swig_destroy__ = _rjmcmc.delete_resultset1dfm
    __del__ = lambda self : None;
    def proposed(self, *args): return _rjmcmc.resultset1dfm_proposed(self, *args)
    def acceptance(self, *args): return _rjmcmc.resultset1dfm_acceptance(self, *args)
    def partitions(self, *args): return _rjmcmc.resultset1dfm_partitions(self, *args)
    def partition_histogram(self, *args): return _rjmcmc.resultset1dfm_partition_histogram(self, *args)
    def partition_location_histogram(self, *args): return _rjmcmc.resultset1dfm_partition_location_histogram(self, *args)
    def x(self, *args): return _rjmcmc.resultset1dfm_x(self, *args)
    def mean(self, li=0): return _rjmcmc.resultset1dfm_mean(self, li)
    def median(self, li=0): return _rjmcmc.resultset1dfm_median(self, li)
    def mode(self, li=0): return _rjmcmc.resultset1dfm_mode(self, li)
    def credible_min(self, li=0): return _rjmcmc.resultset1dfm_credible_min(self, li)
    def credible_max(self, li=0): return _rjmcmc.resultset1dfm_credible_max(self, li)
    def global_parameter(self, *args): return _rjmcmc.resultset1dfm_global_parameter(self, *args)
    def misfit(self, *args): return _rjmcmc.resultset1dfm_misfit(self, *args)
    def __init__(self): 
        """__init__(resultset1dfm self) -> resultset1dfm"""
        this = _rjmcmc.new_resultset1dfm()
        try: self.this.append(this)
        except: self.this = this
resultset1dfm_swigregister = _rjmcmc.resultset1dfm_swigregister
resultset1dfm_swigregister(resultset1dfm)


def regression_single1d(*args):
  """
    regression_single1d(dataset1d dataset, int burnin=10000, int total=50000, int max_order=5, int xsamples=100, 
        int ysamples=100, double credible_interval=0.95) -> resultset1d
    """
  return _rjmcmc.regression_single1d(*args)

def regression_single1d_sampled(*args):
  """
    regression_single1d_sampled(dataset1d dataset, PyObject * callback, int burnin=10000, int total=50000, int max_order=5, 
        int xsamples=100, int ysamples=100, double credible_interval=0.95) -> resultset1d
    """
  return _rjmcmc.regression_single1d_sampled(*args)

def regression_part1d_zero(*args):
  """
    regression_part1d_zero(dataset1d dataset, double pd, int burnin=10000, int total=50000, int max_partitions=20, 
        int xsamples=100, int ysamples=100, double credible_interval=0.95) -> resultset1d
    """
  return _rjmcmc.regression_part1d_zero(*args)

def regression_part1d_natural(*args):
  """
    regression_part1d_natural(dataset1d dataset, double pv, double pd, int burnin=10000, int total=50000, int max_partitions=20, 
        int xsamples=100, int ysamples=100, double credible_interval=0.95) -> resultset1d
    """
  return _rjmcmc.regression_part1d_natural(*args)

def regression_part1d(*args):
  """
    regression_part1d(dataset1d dataset, double pd, int burnin=10000, int total=50000, int max_partitions=20, 
        int max_order=5, int xsamples=100, int ysamples=100, double credible_interval=0.95) -> resultset1d
    """
  return _rjmcmc.regression_part1d(*args)

def regression_part1d_sampled(*args):
  """
    regression_part1d_sampled(dataset1d dataset, PyObject * callback, double pd, int burnin=10000, int total=50000, 
        int max_partitions=20, int max_order=5, int xsamples=100, int ysamples=100, 
        double credible_interval=0.95) -> resultset1d
    """
  return _rjmcmc.regression_part1d_sampled(*args)

def forwardmodel_part1d(*args):
  """
    forwardmodel_part1d(PyObject * local_parameters, PyObject * global_parameters, PyObject * loglikelihood_cb, 
        double minx, double maxx, double pd, int burnin=10000, int total=50000, 
        int max_partitions=20, int xsamples=100, int ysamples=100, double credible_interval=0.95) -> resultset1dfm
    """
  return _rjmcmc.forwardmodel_part1d(*args)
# This file is compatible with both classic and new-style classes.


