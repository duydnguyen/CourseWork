# """
# Debugging demo using pdb
# """

x = "three"
y= -22

def f(z):
    x = z+10
#    import pdb; pdb.set_trace()
    return x

y = f(x)

print "x = ", x
print "y = ", y
