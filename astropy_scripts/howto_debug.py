#Use:
# python howto_debug.py

a = 42
print(a)
b = 2 * a
import IPython; IPython.embed() #before error line
#code is executed until this position,
#then I can check interactively the variables, right before the error occurs, and interactively play to find a solution to the error
#%whos #list of variables
#type(a)
#a.__dict__
1/0 #this gives an error
print(b)
