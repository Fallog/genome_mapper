import y_dc3 as y
import dc3 as s

seq = "abcabcacab"
y.dc3(seq)

# Simon

int_seq = s.add_sentinel_numbers(s.str_to_ascii(seq))
s.dc3(int_seq)
