#test importing script
import sys
sys.path.insert(0,'./../../ImageTools')
import imageTools as IT

print 'test'
test=IT.calculateLineLength([1,1],[2,2])