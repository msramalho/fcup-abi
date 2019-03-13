#!/usr/bin/env python3
# -*- coding: utf-8 -*-

########################### ##################
# Bioinformatics @ DCC-FCUP
# Exercises for class 1
########################### ##################

# Reads a value of temperature in Celsius degrees (oC) and converts it into a temperature in Fahrenheit degrees (oF): F = 32 + 9C/5.
c = int(input("Temperature - Celsius: "))
f = 32 + 9 * c / 5.0
print("Temperature in Fahrenheit:", f)

#Reads a string and converts it to capital letters, printing the result.
seq = input("Insert string: ")
seq1 = seq.upper()
print("String in capital letters:", seq1)

#Adapt the previous program to read the string from a file, whose name is entered by the user.

file_name = input("File name:")
with open (file_name) as fh:
    for line in fh:
        print (line.upper())
fh.close()

#Reads a string and check if it is a palindrome, i.e. if it reads the same when it is reversed. Implement different versions using functions over strings, and cycles (for/while).
str_test = input("String to check:")
str_test_rev = str_test[::-1]
if str_test == str_test_rev:
    print ("string " + str_test + " is palindrome")
else:
    print ("string is not palindrome")

#Reads several positive numerical values from the standard input (until a negative value is found), and calculates the largest and the smallest value.
lst_values = []
val = int(input("insert value:"))
while (val > 0):
    lst_values.append(val)
    val = int(input("insert value:"))
lst_values.sort()
min_val = lst_values[1]
max_val = lst_values[-1]
print ("min: " + str(min_val))
print ("max: " + str(max_val))

#Read an alphanumeric string (long enough) from the input and report the frequencies of each symbol in the sequence.
long_seq = input("insert long sequence:")
freq_hash = {}
for c in long_seq:
    if c in freq_hash:
        freq_hash[c] += 1
    else:
        freq_hash[c] = 1

for k in freq_hash.keys():
    print (str(k) + " "+ str(round(freq_hash[k]/len(long_seq)*100,2)))
 
# Given the two lower sides of a right-angled triangle calculate its hypotenuse.
import math
a = int(input("Length of side A: "))
b = int(input("Length of side B: "))
c = math.sqrt(a ** 2 + b ** 2)
print("Length of hypotenuse:", c)

#Given a numeric interval defined by a lower and upper bound, calculate the sum of all integers included in that interval.
a = int(input("Lower limit of interval: "))
b = int(input("Upper limit of interval: "))
s = 0
for i in range(a, b + 1):
    s += i

print("Sum of integer values in the interval", s)

# defined as a function
def sumInterval(low, upp):
    s = 0
    for i in range(low, upp + 1):
        s += i
    return (s)
