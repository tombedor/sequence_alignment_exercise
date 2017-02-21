import random
letters = ['A', 'C', 'T', 'G']

string = ""
length = 1000
for i in range(length):
    string += letters[random.randrange(0,4)]

print "Answer:"
print string

for i in range(10):
    print ">Frag" + str(i)
    start_i = i*length/10
    end_i = start_i + length / 5

    print string[start_i:end_i]
    


