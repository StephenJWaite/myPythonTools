#Script is designed for testing spline functions
import numpy as np
import imageTools as IT
import matplotlib.pyplot as plt

import numpy as np
print 'running test script'
#xn = [0,1,2,3,4,5,6,7]
#yn = [0,0,0.1,1,2.5,3.2,1,0]
xn = [0,1,2,3]
yn = [0,1,3,-1]
n=np.size(xn)
print 'Test paramaters:'
print '\t-n',n
print '\t-xn',xn
print '\t-yn',yn

[a,b,c,d]=IT.basicCubicSpline(n,xn,yn)

print 'a:',np.shape(a)
print 'b:',np.shape(b)
print 'c:',np.shape(c)
#print 
print 'd:',np.shape(d)

check=IT.basicCubicSplinePlot(a,b.tolist(),c.tolist(),d.tolist(),xn)
if check:
	print 'Test Script Completed, no errors'

plt.plot(xn,yn,'x')
plt.show()

