ncdisp is a simplified "package" (to be packaged and distributed) to be put in the PYTHONPATH to easily load netcdf files and draw 2D maps, lines or vectors of variables.  
    
**Requirements:**  
matplotlib   
basemap  
numpy  
netCDF4  
Python 2.7   
  
If you find any bugs, please do not hesitate to contact the author.  
  
**How to use:    
Example:   
<pre>
nc=ncload('sample')  
u, v=nc.get(['samplevar1','samplevar2'])  
ax,m,im,cb=nc.map(v,loc=['K',1,'O',1],p='nplaea',intp='y',**{'bm':{'boundinglat':20},'map':{'cmap':plt.get_cmap('bwr')},'cb':{'boundaries':range(5)}})  
ax=nc.line([u,v],['I',12,'J',5,'K',1],plotargs=['ro','b-'],**{'ln':{'lw':3},'plot':[{},{'maker':'^'}]'leg':{'loc':1}})
ax,m,Q,qk=nc.vec(u,v,['K',1],p='nplaea',**{'bm':{'boundinglat':20},'q':{'color':'r'},'qk':{'Y':-0.1}})
</pre>   
