#!/usr/bin/env python
from PythonTools import *
import netCDF4

class ncload(ncdisp):
  def __init__(self,fname):
    """
    Load netcdf file.
  
    This function will load a netcdf file using Netcdf4 standard (absolute address or relative address supported). 
    
    Instance attributes:
    -------  
    self.nc             : current netcdf file  
    self.currentdims    : dimensions of the current netcdf file  
    self.currentdimkeys : dimension names of the current netcdf file  
    """
    if fname[-3:]=='.nc':
      ncfile=fname
    else:
      ncfile=fname+'.nc'
#    if True:
    try:
      self.nc=Dataset(ncfile,'r')
      self.name=ncfile
      self.dims,self.dimkeys=self._arrangedims()
#    else:
    except:
      print '>>>>> '+ncfile+' not exist!'
  
  def sh(self,varname=None):
    """
    Show the dimensions of nc file.
  
    Attention:
    The listed dimensions are not garanteed in the right order! Please check the shape of variables before really manipulating.
    Parameters
    ----------
    varname : string, default: None
           If varname is specified, the dimensions of the nc file and the variable are listed.
           If varname is not given, the dimensions of the nc file and all the variables will be list.
    """
    nc=self.nc
    dimkeys=self.dimkeys
    if True:
      varkeys=nc.variables.keys()
      dimkeys=[dimkeys[chr(i+73)] for i in range(len(dimkeys))]
      print dimkeys
      longnames=[np.min([30,len(nc.variables[key].long_name)]) for key in varkeys if key not in dimkeys if 'long_name' in nc.variables[key].ncattrs()]
      if len(longnames)==0:longnames=[0]
      print '>>>>>'
      toprint=['%-20s'%'Dimensions:   ','title   '+' '*(np.max(longnames))]+[chr(k+73)+'     ' for k in range(len(dimkeys))]
      print ''.join(toprint)
      for i,key in enumerate(dimkeys):
        print '  %-18.18s'%(key+':  '),' '*(np.max(longnames)+6),'...   '*i+str('%-6i'%len(nc.dimensions[key]))+'...   '*(len(dimkeys)-i-1)
      print '%-20s'%'Variables:   '+'*'*(6*len(dimkeys)+8+np.max(longnames))
      if not varname:
        for i,key in enumerate(varkeys):
          vardim=nc.variables[key].dimensions
          s=['.','.','.',' ',' ',' ']*len(dimkeys)
          for k,dim in enumerate(vardim):
            index=dimkeys.index(dim)
            s[6*index:6*index+6]=list(str('%-6i'%(nc.variables[key].shape[k])))
          longnameformat='%-'+str(np.max(longnames))+'.'+str(np.max(longnames))+'s'
          if longnames==[0]:
            print '  %-18.18s'%(key+':   '),' '*6,''.join(s)
          else:
            if 'long_name' in nc.variables[key].ncattrs():
              print '  %-18.18s'%(key+':   '),longnameformat%nc.variables[key].long_name+' '*6,''.join(s)
            else:
              print '  %-18.18s'%(key+':   '),' '*(np.max(longnames)+6),''.join(s)
        print '>>>>>'
      else:
        vardim=nc.variables[varname].dimensions                 
        s=['.','.','.',' ',' ',' ']*len(dimkeys)            
        for k,dim in enumerate(vardim):                     
          index=dimkeys.index(dim)                          
          s[6*index:6*index+6]=list(str('%-6i'%(nc.variables[varname].shape[k])))
        longnameformat='%-'+str(np.max(longnames))+'.'+str(np.max(longnames))+'s'
        if longnames==[0]:                                  
          print '  %-18.18s'%(varname+':   '),' '*6,''.join(s)  
        else:                                               
          if 'long_name' in nc.variables[varname].ncattrs():    
            print '  %-18.18s'%(varname+':   '),longnameformat%nc.variables[varname].long_name+' '*6,''.join(s)
          else:                                             
            print '  %-18.18s'%(varname+':   '),' '*(np.max(longnames)+6),''.join(s)
        print '>>>>>' 
  
  def _getvar(self,vname):
    """
    Choose a variable in the nc file.
    
    This function will read the variable named 'varname'.
   
    Returns
    -------
    var : netCDF.Variable object
    """
    return self.nc.variables[vname]

  def get(self,*vnames):
    """
    Choose a variable or a list of variales.
    
    Parameters
    ----------
    varnames : list or string
           Names of variables. Could be string (one variable) or list (several variables)

    Returns
    -------
    varlist: netCDF4.Variable objects
    """
    varlist=[self._getvar(vname) for vname in vnames]
    self._refvar=varlist[0]
    return varlist[0] if len(vnames)==1 else varlist
  
  def _arrangedims(self):
    """
    Arrange the dimensions in the nc file.
    
    Return:
    -------  
    dims : keys from 'I' to map the dimensions of the nc file
    dimkeys : keys from 'I' to map the dimension names of the nc file
    """
    nc=self.nc
    dims={}
    dimkeys={}
    ncdimkeys=nc.dimensions.keys()
    for key in ncdimkeys:
      if 'longitude' in key:
        lonindex=ncdimkeys.index(key)
        ncdimkeys.insert(0,ncdimkeys.pop(lonindex))
        continue
      if ('longitude' not in key) and ('lon' in key):
        lonindex=ncdimkeys.index(key)
        ncdimkeys.insert(0,ncdimkeys.pop(lonindex))
        continue
      if 'latitude' in key:
        latindex=ncdimkeys.index(key)
        ncdimkeys.insert(0,ncdimkeys.pop(latindex))
        continue
      if ('latitude' not in key) and ('lat' in key):
        latindex=ncdimkeys.index(key)
        ncdimkeys.insert(0,ncdimkeys.pop(latindex))
        continue
    for i,key in enumerate(ncdimkeys):
      if key in nc.variables.keys():
        dims[chr(73+i)]=nc.variables[key][:] 
      else:dims[chr(73+i)]=np.arange(1,len(nc.dimensions[key])+1)
      dimkeys[chr(73+i)]=key
    return dims,dimkeys
  
  def close(self):
    self.nc.close()
