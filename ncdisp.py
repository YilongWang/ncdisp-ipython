#!/usr/bin/env python
import numpy as np
import copy
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import inspect
import matplotlib.ticker as ticker 
from netCDF4 import Dataset

class ncdisp:      
  def map(self,var,y1=None,y2=None,x1=None,x2=None,y='I',x='J',loc=[],\
            m=None,p='cyl',edg='n',face='y',intp='n',cbar='y',ref=None,**kwargs):
    """
    Display the geographic map.
  
    Parameters
    ----------
    var : variable to display;
          Could be a extracted matrix or a string indicating the variable name
    y, x : [ 'I' | 'J' | ... | 'i' | 'j' | ... |], default: 'I', 'J'
           Dimension of the data. 
           If lowercase, 'y1' ('x1') and 'y2' ('x2') should be position, e.g. latitude/longitude/time
           Otherwise, 'y1' ('x1') and 'y2' ('x2') should be indexes
    y1, y2, x1, x2 : map boundaries, default: None
    loc : [['I',1...] | ['i',1...]], default: []
           Locator of the variable on other dimensions. 
           When the letters are uppercases, the number after it should be pixel index;
           When the letters are lowercases, the number after it should be exact position (e.g. latitude/longitude/time)
    m : basemap instance, default: None
           If m is given, draw the map on m; otherwise, draw on a new basemap.
  
    p (short for projection) : [ None | projection for mpl_toolkits.basemap.Basemap], default: 'cyl'
           If 'none', common pcolormesh or imshow will be invoked 
    edge : [ 'y' | 'n' ], default: 'n'
           Whether (y) of not (n) to draw contour lines
    face : [ 'y' | 'n' ], default: 'y'
           Whether (y) of not (n) to draw contour surface
    cbar : [ 'y' | 'n' ], default: 'y'
           If 'y', draw the colormap attached to the map.
    ref : netCDF4.Variable object, default: None
           Variable for reference of dimensions. If None, self._refvar will be used.
    **kwargs : 
           Dict with up to four inferior dict (optional) are available:
              'bm' : **bmkwargs, passed on to mpl_toolkits.basemap.Basemap
              'map' : **mapkwargs, passed on to matplotlib.pyplot.contour(f) or matplotlib.pyplot.pcolormesh or matplotlib.pyplot.imshow
              'cb' : **cbkwargs, passed on to matplotlib.pyplot.colorbar

    Returns  
    -------  
    ax : axes  
    m  : object of basemap or axes  
    im : image  
    cb : colorbar  

    Example:
    nc=ncload('sample')
    v=nc.get('samplevar')
    ax,m,im,cb=nc.map(v,loc=['K',1,'O',1],p='nplaea',intp='y',bm={'boundinglat':20},map={'cmap':plt.get_cmap('bwr')},cb={'boundaries':range(5)}})
    """
    nc=self.nc
    dims=self.dims
    dimkeys=self.dimkeys
    
    bmspec=inspect.getargspec(self._defaultbmargs)
    bmkwargs=('bm' in kwargs and kwargs['bm']) or {}
    mapspec=inspect.getargspec(self._defaultmapargs)
    mapkwargs=('map' in kwargs and kwargs['map']) or {}
    cbkwargs=('cb' in kwargs and kwargs['cb']) or {}
    dic=dict(zip(bmspec[0],bmspec[-1]))
    for key in dic.keys():
      if key not in bmkwargs:bmkwargs[key]=dic[key]
    dic=dict(zip(mapspec[0],mapspec[-1]))
    for key in dic.keys():
      if key not in mapkwargs:mapkwargs[key]=dic[key]

    if loc:
      locdims=[loc[i*2] for i in range(len(loc)/2)]
      locposs=[loc[i*2+1] for i in range(len(loc)/2)]
      index=[slice(0,k) for k in var.shape]
      if isinstance(var,netCDF4.Variable):ref=var
      elif ref==None:ref=self._refvar
      for d in locdims:index[ref.dimensions.index(dimkeys[d.upper()])]=locposs[locdims.index(d)]
      var=var[index]
    
    coordy=dims[y.upper()]
    coordx=dims[x.upper()]
    if coordy[0]>coordy[-1]:
      var=var[::-1]
      coordy=coordy[::-1]
      if y.isupper():
        ytemp=y1
        if y2!=None:y1=len(coordy)-y2
        else:y1=None
        if ytemp!=None:y2=len(coordy)-ytemp
        else:y2=None
    if coordx[0]>coordx[-1]:
      var=var[:,::-1]
      coordx=coordx[::-1] 
      if x.isupper():
        xtemp=x1
        if x2!=None:x1=len(coordx)-x2
        else:x1=None
        if xtemp!=None:x2=len(coordx)-xtemp
        else:x2=None
  
    if 'lat' in dimkeys[y.upper()].lower():
      latmin=np.max([-90.,coordy[0]-(coordy[1]-coordy[0])/2.])
      latmax=np.min([90.,coordy[-1]+(coordy[-1]-coordy[-2])/2.])
    elif 'lon' in dimkeys[y.upper()].lower():
      coordy=np.hstack((coordy-360.,coordy,coordy+360.))
      coordymask=(coordy>=-180.)*(coordy<=180.)
      coordy=coordy[coordymask]
      coordy=coordy[(coordy[0]==coordy[1] and 1) or 0:(coordy[-1]==coordy[-2] and -1) or None]
      var=np.ma.vstack((var,var,var))
      var=var[coordymask,:]
      var=var[(coordy[0]==coordy[1] and 1) or 0:(coordy[-1]==coordy[-2] and -1) or None,:]
      latmin=np.max([-180.,coordy[0]-(coordy[1]-coordy[0])/2.])
      latmax=np.min([180.,coordy[-1]+(coordy[-1]-coordy[-2])/2.])
    else:
      latmin=latmax=np.nan
    if 'lat' in dimkeys[x.upper()].lower():
      lonmin=np.max([-90.,coordx[0]-(coordx[1]-coordx[0])/2.])
      lonmax=np.min([90.,coordx[-1]+(coordx[-1]-coordx[-2])/2.])
    elif 'lon' in dimkeys[x.upper()].lower():
      coordx=np.hstack((coordx-360.,coordx,coordx+360.))
      coordxmask=(coordx>=-180.)*(coordx<=180.)
      coordx=coordx[coordxmask]
      coordx=coordx[(coordx[0]==coordx[1] and 1) or 0:(coordx[-1]==coordx[-2] and -1) or None]
      var=np.ma.hstack((var,var,var))
      var=var[:,coordxmask]
      var=var[:,(coordx[0]==coordx[1] and 1) or 0:(coordx[-1]==coordx[-2] and -1) or None]
      lonmin=np.max([-180.,coordx[0]-(coordx[1]-coordx[0])/2.])
      lonmax=np.min([180.,coordx[-1]+(coordx[-1]-coordx[-2])/2.])
    else:
      lonmin=lonmax=np.nan

    if p in ['nplaea','ortho']:
      crnryy=(coordy[:-1]+coordy[1:])/2.
      crnryy=np.insert(crnryy,0,coordy[0]-(coordy[1]-coordy[0])/2.)
      crnryy=np.append(crnryy,coordy[-1]+(coordy[-1]-coordy[-2])/2.)
      crnrxx=(coordx[:-1]+coordx[1:])/2.
      crnrxx=np.insert(crnrxx,0,coordx[0]-(coordx[1]-coordx[0])/2.)
      crnrxx=np.append(crnrxx,coordx[-1]+(coordx[-1]-coordx[-2])/2.)
    else:
      crnryy=(coordy[:-1]+coordy[1:])/2.
      crnryy=np.insert(crnryy,0,np.nanmax([latmin,coordy[0]-(coordy[1]-coordy[0])/2.]))
      crnryy=np.append(crnryy,np.nanmin([latmax,coordy[-1]+(coordy[-1]-coordy[-2])/2.]))
      crnrxx=(coordx[:-1]+coordx[1:])/2.
      crnrxx=np.insert(crnrxx,0,np.nanmax([lonmin,coordx[0]-(coordx[1]-coordx[0])/2.]))
      crnrxx=np.append(crnrxx,np.nanmin([lonmax,coordx[-1]+(coordx[-1]-coordx[-2])/2.]))
  
    if y.isupper():
      if y1==None:y1=0
      if y2==None:y2=len(dims[y.upper()])-1
      py=[y1,y2]
    else:
      if y1==None:y1=coordy[0]
      if y2==None:y2=coordy[-1]
      py=np.searchsorted(crnryy,[y1,y2])-1
    if x.isupper():
      if x1==None:x1=0
      if x2==None:x2=len(dims[x.upper()])-1
      px=[x1,x2]
    else:
      if x1==None:x1=coordx[0]
      if x2==None:x2=coordx[-1]
      px=np.searchsorted(crnrxx,[x1,x2])-1

    mapvar=var[py[0]:py[-1]+1,px[0]:px[-1]+1]
    if intp=='n':
      mapboundy=[crnryy[py[0]],crnryy[py[-1]+1]]
      mapboundx=[crnrxx[px[0]],crnrxx[px[-1]+1]]
    else:
      mapboundy=[coordy[py[0]],coordy[py[-1]]]
      mapboundx=[coordx[px[0]],coordx[px[-1]]]
  
    ax=plt.gca()
    if not m:
      if p:
        if p=='cyl':
          if not set(['llcrnrlat','llcrnrlon','urcrnrlat','urcrnrlon']).intersection(bmkwargs.keys()):
            args=[mapboundx[0],mapboundy[0],mapboundx[1],mapboundy[1]]
            m=Basemap(*args,**bmkwargs)                                
          else:
            m=Basemap(**bmkwargs)
        else:
          m=Basemap(projection=p,**bmkwargs)
  
  #      m.drawcoastlines()
  #      ax.axis('off')
        if p in ['nplaea','ortho']:
          m.drawparallels(np.arange(-80.,81.,40.))
          m.drawmeridians(np.arange(-180.,181.,60.),labels=[1,1,1,1])
        else:
          m.drawparallels(np.arange(-80.,81.,40.),labels=[1,0,0,0])
          m.drawmeridians(np.arange(-180.,181.,60.),labels=[0,0,0,1])
      else:
        m=ax

    if p:
      if intp=='n':
        X,Y=np.meshgrid(crnrxx[px[0]:px[-1]+2],crnryy[py[0]:py[-1]+2])
        dispx,dispy=m(X,Y)
        im=m.pcolormesh(dispx,dispy,mapvar,**mapkwargs)
      else:
        X,Y=np.meshgrid(coordx[px[0]:px[-1]+1],coordy[py[0]:py[-1]+1])
        dispx,dispy=m(X,Y)
        if 'ticks' in cbkwargs:cmin=cbkwargs['ticks'][0];cmax=cbkwargs['ticks'][-1]
        else:
          if 'vmin' in mapkwargs:cmin=mapkwargs['vmin']
          else:cmin=np.nanmin(mapvar)
          if 'vmax' in mapkwargs:cmax=mapkwargs['vmax']
          else:cmax=np.nanmax(mapvar)
        locator = ticker.MaxNLocator(256)
        locator.create_dummy_axis()
        locator.set_bounds(cmin, cmax) 
        levs = locator()
        if face=='y':   
          if 'colorsf' in mapkwargs:
            im=m.contourf(dispx,dispy,mapvar,levs,colors=mapkwargs['colorsf'],**mapkwargs)
          else:      
            im=m.contourf(dispx,dispy,mapvar,levs,**mapkwargs)
        if edg=='y':
          if 'colorsl' in mapkwargs and 'cmap' in mapkwargs:
            im=m.contour(dispx,dispy,mapvar,colors=mapkwargs['colorsl'],**mapkwargs)
          else:
            im=m.contour(dispx,dispy,mapvar,**mapkwargs)
    else:          # if not p
      if intp=='n':
        im=m.pcolormesh(crnrxx[px[0]:px[-1]+2],crnryy[py[0]:py[-1]+2],mapvar,**mapkwargs)
        m.set_xlim(crnrxx[px[0]],crnrxx[px[-1]+1])
        m.set_ylim(crnryy[py[0]],crnryy[py[-1]+1])
      else:
        X,Y=np.meshgrid(coordx[px[0]:px[-1]+1],coordy[py[0]:py[-1]+1])
        if 'ticks' in cbkwargs:cmin=cbkwargs['ticks'][0];cmax=cbkwargs['ticks'][-1]
        else:
          if 'vmin' in mapkwargs:cmin=mapkwargs['vmin']
          else:cmin=np.nanmin(mapvar)
          if 'vmax' in mapkwargs:cmax=mapkwargs['vmax']
          else:cmax=np.nanmax(mapvar)
        locator = ticker.MaxNLocator(256)
        locator.create_dummy_axis()
        locator.set_bounds(cmin, cmax) 
        levs = locator()
        if face=='y':   
          if 'colorsf' in mapkwargs:
            im=m.contourf(dispx,dispy,mapvar,levs,colors=mapkwargs['colorsf'],**mapkwargs)
          else:      
            im=m.contourf(dispx,dispy,mapvar,levs,**mapkwargs)
        if edg=='y':
          if 'colorsl' in mapkwargs and 'cmap' in mapkwargs:
            im=m.contour(dispx,dispy,mapvar,colors=mapkwargs['colorsl'],**mapkwargs)
          else:
            im=m.contour(dispx,dispy,mapvar,**mapkwargs)
  
    if cbar=='y':
      cbkwargs_tmp={}
      for key in cbkwargs:
        if key!='ticklabels':cbkwargs_tmp[key]=cbkwargs[key]
      cb=plt.colorbar(im,**cbkwargs_tmp)
      if 'ticklabels' in cbkwargs:
        if 'orientation' in cbkwargs and cbkwargs['orientation']=='horizontal':
          cb.ax.set_xticklabels(cbkwargs['ticklabels'])
        else:
          cb.ax.set_yticklabels(cbkwargs['ticklabels'])
    else:cb=0
    
    return ax,m,im,cb
  
  def _defaultbmargs(lat_0=90,lon_0=0,boundinglat=1):
    pass

  def _defaultmapargs(zorder=0):
    pass
  
  def line(self,vs=[],loc=[],orientation='h',invertyaxis='n',legend='on',refs=[],plotargs=[],**kwargs):
    """
    Draw the line.
    Attention: 
           Global variables 'ax', 'lines', 'legends' will be overwritten or set up.
    Parameters
    ----------
    vs : [...], default: None
           Names for variables to plot. Must be a list of str obejects.
    loc : [['I',1...] | ['i',1...]], default: None
           Position of the varibles. When the letters are uppercases, the number after it should be pixel index;
           When the letters are lowercases, the number after it should be position (e.g. latitude/longitude/time)
    orientation : [ 'h' | 'v'], default: 'h'
           Orientation of the lines.
    invertyaxis : [ 'n' | 'y' ], default: 'n'
           Choose whether to invert the y axis.
    legend : [ 'on' | 'off' ], default: 'on'
           If 'on', add the legend; else, there is no legend
    ref : netCDF4.Variable objects, default: None
           Variables for reference of dimensions.
           If None, the first variable activated by last call of nc.get() will be used.
    *plotargs :
           List of *args that will be passed on to matplotlib.pyplot
    **kwargs : 
           Dict with up to three inferior dict (optional) are available:
              'ln'  : **linekwargs, passed on to matplolib.pyplot.plot  
                       Common features  
              'plot' : a list of **plotkwargs, passed on to matplolib.pyplot.plot  
                       Specific features  
              'leg' : **legendkwargs, passed on to matplolib.pyplot.legend  
    
    Returns  
    -------  
    ax : axes  

    Example:
    nc=ncload('sample')
    v1,v2=nc.get(['samplevar1','samplevar2'])
    v3=v1[:]+v2[:]
    nc.line([v1,v2,v3],['I',12,'J',5],'y',ref=[v1,v2,v1],plotargs=['ro','b-','y'],ln={'lw':3},plot=[{},{},{'maker':'^'}],leg={'loc':1})
    """
    nc=self.nc
    dims=self.dims
    dimkeys=self.dimkeys

    lineskwargs=('ln' in kwargs and kwargs['ln']) or {}
    plotkwargs=('plot' in kwargs and kwargs['plot']) or []
    legendkwargs=('legend' in kwargs and kwargs['legend']) or {}

    locdims=[loc[i*2] for i in range(len(loc)/2)]
    locposs=[loc[i*2+1] for i in range(len(loc)/2)]
    
    ax=plt.gca()

    ylim1=0.
    ylim2=0.
    for v in vs:
      if isinstance(v,netCDF4.Variable):ref=v
      elif refs==[]:ref=self._refvar
      elif len(refs)==1:ref=refs[0]
      elif vs.index(v)>=len(ref):ref=refs[-1]
      else:ref=refs[vs.index(v)]
      index=[range(k) for k in v.shape]
      axis=[]
      invert='n'
      for d in locdims:
        if type(locposs[locdims.index(d)])==list or type(locposs[locdims.index(d)])==np.ndarray:
#          if str(dims[d.upper()].dtype)[1]!='S':axis=dims[d.upper()][oneindex]
#          else:axis=lineposs[linedims.index(d)]
          if d.isupper():
            start,end=locposs[locdims.index(d)]
            if dims[d.upper()][-1]<dims[d.upper()][0] and d.upper()!='I' and d.upper()!='J':invert='y'
          else:
            if dims[d.upper()][-1]>dims[d.upper()][0]:
              start,end=np.searchsorted(dims[d.upper()],locposs[locdims.index(d)])
            else:
              end,start=len(dims[d.upper()])-np.searchsorted(dims[d.upper()][::-1],locposs[locdims.index(d)])
              if d.upper()!='I' and d.upper()!='J':invert='y'
          oneindex=slice(np.max([0,start-1]),np.min([len(dims[dim.upper()]),end+1]))
          index[ref.dimensions.index(dimkeys[d.upper()])]=oneindex
          axis=dims[d.upper()][index]
        else:
          if d.isupper():                                       
            oneindex=locposs[locdims.index(d)]
            if dims[d.upper()][-1]<dims[d.upper()][0] and d.upper()!='I' and d.upper()!='J':invert='y'
          else:                      
            if dims[d.upper()][-1]>dims[d.upper()][0]:
              oneindex=np.searchsorted(dims[d.upper()],locposs[locdims.index(d)])
            else:                    
              oneindex=len(dims[d.upper()])-np.searchsorted(dims[d.upper()][::-1],locposs[locdims.index(d)])
              if d.upper()!='I' and d.upper()!='J':invert='y'
          index[ref.dimensions.index(dimkeys[d.upper()])]=oneindex
      if axis==[]:    # that means the dim not given is used to draw the line
        for dim in ref.dimensions:
          for key in dimkeys.keys():
            if dim==dimkeys[key]:break
          if (key not in locdims and v.shape[ref.dimensions.index(dim)]!=1) and\
             (key.lower() not in locdims and v.shape[ref.dimensions.index(dim)]!=1):
            axis=dims[key][:];break
          else:continue
    
      value=v[tuple(index)]
      if axis==[]: axis=range(1,len(value)+1)
      if orientation=='h':x=axis;y=value
      else: x=value;y=axis
      ylim1=np.min([np.min(y),ylim1])
      ylim2=np.max([np.max(y),ylim2])
      label=ref._name+' '
      for i in range(len(locdims)):
        if type(locposs[i])==list:
          label=label+locdims[i]+'='+str(locposs[i][0])+'-'+str(locposs[i][-1])+' '
        else:
          label=label+locdims[i]+'='+str(locposs[i])+' '
      if plotargs:
        linekwargs=lineskwargs
        if plotkwargs:linekwargs.update(plotkwargs[vs.index(v)])
        ln=ax.plot(x,y,plotargs[vs.index(v)],label=label,**linekwargs)
      else:
        linekwargs=lineskwargs
        if plotkwargs:linekwargs.update(plotkwargs[vs.index(v)])
        ln=ax.plot(x,y,label=label,**linekwargs) 
    if orientation=='h':
      if invert=='y' or invertyaxis=='y':ax.invert_xaxis()
    else:
      if invert=='y' or invertyaxis=='y':ax.set_ylim(ylim1-(ylim2-ylim1)/20.,ylim2);ax.invert_yaxis()
      else:ax.set_ylim(ylim1,ylim2+(ylim2-ylim1)/20.)
    if legend=='on':
      legends=pylab.legend(**legendkwargs)
  
    return ax
   
  def vec(self,u,v,y1=None,y2=None,x1=None,x2=None,y='I',x='J',loc=[],\
            m=None,ny=20,nx=20,key='y',p='cyl',ref=None,**kwargs):
    """
    Display the vector (e.g., wind field).
  
      Attention: 
           Global variables 'ax', 'm', 'Q', 'qk' will be set up or updated.
  
      *ax*  : axes of the image
      *m*   : Basemap object
      *Q*   : quiver object
      *qk*  : quiverkey object
    
    Parameters
    ----------
    u : variable for horizontal direction
    v : variable for vertical direction
    y, x : [ 'I' | 'J' | ... | 'i' | 'j' | ... |], default: 'I', 'J'
           Dimension of the data. 
           If lowercase, 'y1' ('x1') and 'y2' ('x2') should be position, e.g. latitude/longitude/time
           Otherwise, 'y1' ('x1') and 'y2' ('x2') should be indexes
    y1, y2, x1, x2 : map boundaries, default: None
    loc : [['I',1...] | ['i',1...]], default: []
           Locator of the variable on other dimensions. 
           When the letters are uppercases, the number after it should be pixel index;
           When the letters are lowercases, the number after it should be exact position (e.g. latitude/longitude/time)
    m : basemap instance, default: None
           If m is given, draw the map on m; otherwise, draw on a new basemap.
    ny, nx : int, default: 20, 20
           numbers of rows and columns of vectors
  
    p (short for projection) : [ None | projection for mpl_toolkits.basemap.Basemap], default: 'cyl'
           If 'none', common pcolormesh or imshow will be invoked 
    ref : netCDF4.Variable object, default: None
           Variable for reference of dimensions. If None, self._refvar will be used.
    **kwargs : 
           Dict with up to four inferior dict (optional) are available:
              'bm' : **bmkwargs, passed on to mpl_toolkits.basemap.Basemap
              'q'  : **caxkwargs, passed on to mpl_toolkits.basemap.quiver
              'qk' : **mapkwargs, passed on to matplotlib.pyplot.quiverkey

    Returns  
    -------  
    ax : axes  
    m : object of basemap or axes
    Q : quiver  
    qk (optional) : quiverkey  

    Example:
    nc=ncload('sample')
    u,v=nc.get(['uwind','vwind'])
    ax,m,Q,qk=nc.vec(u,v,['K',1],p='nplaea',bm={'boundinglat':20},q={'color':'r'},qk={'Y':-0.1}})
    """
    nc=self.nc
    dims=self.dims
    dimkeys=self.dimkeys
  
    qkspec=inspect.getargspec(self._defaultqkargs)
    quiverkwargs=('q' in kwargs and kwargs['q']) or {}
    qkkwargs=('qk' in kwargs and kwargs['qk']) or {}
    bmspec=inspect.getargspec(self._defaultbmargs)
    bmkwargs=('bm' in kwargs and kwargs['bm']) or {}
    dic=dict(zip(bmspec[0],bmspec[-1]))
    for bmkey in dic.keys():
      if bmkey not in bmkwargs:bmkwargs[bmkey]=dic[bmkey]
    dic=dict(zip(qkspec[0],qkspec[-1]))
    for qkkey in dic.keys():
      if qkkey not in qkkwargs:qkkwargs[qkkey]=dic[qkkey]

    if loc:
      locdims=[loc[i*2] for i in range(len(loc)/2)]
      locposs=[loc[i*2+1] for i in range(len(loc)/2)]
      index=[slice(0,k) for k in v.shape]
      if isinstance(u,netCDF4.Variable):ref=u
      elif isinstance(v,netCDF4.Variable):ref=v
      elif ref==None:ref=self._refvar
      for d in locdims:index[ref.dimensions.index(dimkeys[d.upper()])]=locposs[locdims.index(d)]
      v=v[index]
      u=u[index]

    coordy=dims[y.upper()]
    coordx=dims[x.upper()]
    if coordy[0]>coordy[-1]:
      v=v[::-1]
      u=u[::-1]
      coordy=coordy[::-1]
      if y!='x' and y!='y':
        ytemp=y1
        if y2!=None:y1=len(coordy)-y2
        else:y1=None
        if ytemp!=None:y2=len(coordy)-ytemp
        else:y2=None
    if coordx[0]>coordx[-1]:
      v=v[:,::-1]
      u=u[:,::-1]
      coordx=coordx[::-1] 
      if x!='x' and x!='y':
        xtemp=x1
        if x2!=None:x1=len(coordx)-x2
        else:x1=None
        if xtemp!=None:x2=len(coordx)-xtemp
        else:x2=None
  
    if 'lat' in dimkeys[y.upper()].lower():
      latmin=np.max([-90.,coordy[0]-(coordy[1]-coordy[0])/2.])
      latmax=np.min([90.,coordy[-1]+(coordy[-1]-coordy[-2])/2.])
    elif 'lon' in dimkeys[y.upper()].lower():
      coordy=np.hstack((coordy-360.,coordy,coordy+360.))
      coordymask=(coordy>=-180.)*(coordy<=180.)
      coordy=coordy[coordymask]
      coordy=coordy[(coordy[0]==coordy[1] and 1) or 0:(coordy[-1]==coordy[-2] and -1) or None]
      v=np.ma.hstack((v,v,v))
      v=v[:,coordxmask]
      v=v[:,(coordx[0]==coordx[1] and 1) or 0:(coordx[-1]==coordx[-2] and -1) or None]
      u=np.ma.hstack((u,u,u))
      u=u[:,coordxmask]
      u=u[:,(coordx[0]==coordx[1] and 1) or 0:(coordx[-1]==coordx[-2] and -1) or None]
      latmin=np.max([-180.,coordy[0]-(coordy[1]-coordy[0])/2.])
      latmax=np.min([180.,coordy[-1]+(coordy[-1]-coordy[-2])/2.])
    else:
      latmin=latmax=np.nan
    if 'lat' in dimkeys[x.upper()].lower():
      lonmin=np.max([-90.,coordx[0]-(coordx[1]-coordx[0])/2.])
      lonmax=np.min([90.,coordx[-1]+(coordx[-1]-coordx[-2])/2.])
    elif 'lon' in dimkeys[x.upper()].lower():
      coordx=np.hstack((coordx-360.,coordx,coordx+360.))
      coordxmask=(coordx>=-180.)*(coordx<=180.)
      coordx=coordx[coordxmask]
      coordx=coordx[(coordx[0]==coordx[1] and 1) or 0:(coordx[-1]==coordx[-2] and -1) or None]
      v=np.ma.hstack((v,v,v))
      v=v[:,coordxmask]
      v=v[:,(coordx[0]==coordx[1] and 1) or 0:(coordx[-1]==coordx[-2] and -1) or None]
      u=np.ma.hstack((u,u,u))
      u=u[:,coordxmask]
      u=u[:,(coordx[0]==coordx[1] and 1) or 0:(coordx[-1]==coordx[-2] and -1) or None]
      lonmin=np.max([-180.,coordx[0]-(coordx[1]-coordx[0])/2.])
      lonmax=np.min([180.,coordx[-1]+(coordx[-1]-coordx[-2])/2.])
    else:
      lonmin=lonmax=np.nan

    if p in ['nplaea','ortho']:
      crnryy=(coordy[:-1]+coordy[1:])/2.
      crnryy=np.insert(crnryy,0,coordy[0]-(coordy[1]-coordy[0])/2.)
      crnryy=np.append(crnryy,coordy[-1]+(coordy[-1]-coordy[-2])/2.)
      crnrxx=(coordx[:-1]+coordx[1:])/2.
      crnrxx=np.insert(crnrxx,0,coordx[0]-(coordx[1]-coordx[0])/2.)
      crnrxx=np.append(crnrxx,coordx[-1]+(coordx[-1]-coordx[-2])/2.)
    else:
      crnryy=(coordy[:-1]+coordy[1:])/2.
      crnryy=np.insert(crnryy,0,np.nanmax([latmin,coordy[0]-(coordy[1]-coordy[0])/2.]))
      crnryy=np.append(crnryy,np.nanmin([latmax,coordy[-1]+(coordy[-1]-coordy[-2])/2.]))
      crnrxx=(coordx[:-1]+coordx[1:])/2.
      crnrxx=np.insert(crnrxx,0,np.nanmax([lonmin,coordx[0]-(coordx[1]-coordx[0])/2.]))
      crnrxx=np.append(crnrxx,np.nanmin([lonmax,coordx[-1]+(coordx[-1]-coordx[-2])/2.]))
  
    if y.isupper():
      if y1==None:y1=0
      if y2==None:y2=len(coordy)-1
      py=[y1,y2]
    else:
      if y1==None:y1=coordy[0]
      if y2==None:y2=coordy[-1]
      py=np.searchsorted(crnryy,[y1,y2])-1
    if x.isupper():
      if x1==None:x1=0
      if x2==None:x2=len(coordx)-1
      px=[x1,x2]
    else:
      if x1==None:x1=coordx[0]
      if x2==None:x2=coordx[-1]
      px=np.searchsorted(crnrxx,[x1,x2])-1

    v=v[py[0]:py[-1]+1,px[0]:px[-1]+1]
    u=u[py[0]:py[-1]+1,px[0]:px[-1]+1]
    mapboundy=[crnryy[py[0]],crnryy[py[-1]+1]]
    mapboundx=[crnrxx[px[0]],crnrxx[px[-1]+1]]

    ax=plt.gca()
    if not m:
      if p:
        if p=='cyl':
          if not set(['llcrnrlat','llcrnrlon','urcrnrlat','urcrnrlon']).intersection(bmkwargs.keys()):
            args=[mapboundx[0],mapboundy[0],mapboundx[1],mapboundy[1]]
            m=Basemap(*args,**bmkwargs)                                
          else:
            m=Basemap(**bmkwargs)
        else:
          m=Basemap(projection=p,**bmkwargs)
  
        if p in ['nplaea','ortho']:
          m.drawparallels(np.arange(-80.,81.,40.))
          m.drawmeridians(np.arange(-180.,181.,60.),labels=[1,1,1,1])
        else:
          m.drawparallels(np.arange(-80.,81.,40.),labels=[1,0,0,0])
          m.drawmeridians(np.arange(-180.,181.,60.),labels=[0,0,1,0])
      else:
        m=ax
  
    intervx=v.shape[1]//nx
    intervy=v.shape[0]//ny
    X,Y=np.meshgrid(coordx[px[0]:px[-1]+1:intervx],coordy[py[0]:py[-1]+1:intervy])
    dispx,dispy=m(X,Y)
    Q = m.quiver(dispx,dispy,u[::intervy,::intervx],v[::intervy,::intervx],**quiverkwargs)
#    uproj,vproj,dispx,dispy = m.transform_vector(u[::intervy,::intervx],v[::intervy,::intervx],\
#                        coordx[px[0]:px[-1]+1:intervx],coordy[py[0]:py[-1]+1:intervy],nx,ny,returnxy=True)
#    Q = m.quiver(dispx,dispy,uproj,vproj,**quiverkwargs)
    if key=='y':
      length=int(np.nanmax(np.sqrt(v**2+u**2)))
      if 'U' not in qkkwargs.keys():qkkwargs['U']=length
      if 'label' not in qkkwargs.keys():
        qkkwargs['label']=str(length)+' '+(('units' in ref.ncattrs() and ref.units) or '')
      qk = plt.quiverkey(Q,**qkkwargs)
      return ax,m,Q,qk 
    else: return ax,m,Q

  def _defaultqkargs(X=0.5,Y=-0.1):
    pass
