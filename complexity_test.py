import y_dc3 as y
import dc3 as s
import cProfile

seq = "ACCGTAGGCCATT" * 5000

cProfile.run("y.dc3(seq)")
# y.dc3(seq)

# Simon

int_seq = s.add_sentinel_numbers(s.str_to_ascii(seq))
# s.dc3(int_seq)
cProfile.run("s.dc3(int_seq)")
