import y_dc3 as y
import dc3 as s
import cProfile
import random

# seq = "actgtatgtatatatgatatagaataatagtatat" * 1000
seq = ""
for i in range(3330000):
    seq += random.choice(["a", "c", "g", "t"])

# print(seq)
# print(y.dc3(seq))

cProfile.run("""y.dc3(seq, sorting_algorithm = "stable")""")
# # print(y.dc3(seq, True))
# print(y.dc3(seq))

## 44 sec : 600 000
## 5 minutes : 3M3
