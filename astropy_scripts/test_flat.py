import numpy as np

x = np.arange(1, 7).reshape(2, 3)

print "x"
print x

print ""

print "x.flat[3]"
print x.flat[3]

print ""

print "x.flat"
print x.flat


fl = x.flat

print ""

print "fl"
print fl

print ""

print "type(fl)"
print type(fl)

print ""

print "iterate over fl"
for item in fl:
    print item
