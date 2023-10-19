import subprocess as sp


p1 = sp.Popen(["lsusb"], stdout=sp.PIPE)
p2 = sp.Popen(["grep", "Oxford"], stdin=p1.stdout, stdout=sp.PIPE,)
p1.stdout.close()  # Allow p1 to receive a SIGPIPE if p2 exits.
output = p2.communicate()
print("output")