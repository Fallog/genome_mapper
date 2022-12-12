import numpy as np
import tools


class Sais:
    s_type = True
    l_type = False

    def build_suffix_array(self, txt: np.ndarray, alpha_size=5):
        length = txt.size
        types = self.get_suffix_type(txt, length)
        size_dict = self.get_bucket_sizes(txt, length, alpha_size)
        suffix_array = np.full(length, -1)

        self.poisition_lms_char(txt, length, suffix_array, types, size_dict)
        # print(suffix_array, "1")
        self.induction_sort_l(txt, suffix_array, types, size_dict)
        # print(suffix_array, "2")
        self.induction_sort_s(txt, suffix_array, types, size_dict)
        # print(suffix_array, "SA")

        reduced_txt, vals, current_names = self.reduce(txt, length, suffix_array, types)
        # print(reduced_txt, "y")
        # print(vals, "vals")
        last_suffixe_array = self.build_summary_suffix_array(
            reduced_txt, suffix_array, current_names
        )
        if np.array_equal(last_suffixe_array, suffix_array):
            return last_suffixe_array
        # print(last_suffixe_array, "last LA")
        suffix_array = np.full(suffix_array.size, -1, dtype=int)
        vals = np.take_along_axis(vals, last_suffixe_array, axis=0)
        # print(vals, "after sort")
        self.last_poisition_lms_char(
            txt, length, suffix_array, types, size_dict, last_suffixe_array, vals
        )
        # print(suffix_array)
        self.induction_sort_l(txt, suffix_array, types, size_dict)

        self.induction_sort_s(txt, suffix_array, types, size_dict)
        # print("end")
        return suffix_array

    def last_poisition_lms_char(
        self,
        txt: np.ndarray,
        length: int,
        suffix_array: np.ndarray,
        types_list: np.ndarray,
        bucket_sizes: dict[int],
        last_suffix_array: np.ndarray,
        vals: np.ndarray,
    ):
        suffix_array[0] = length
        bucket_tails = self.get_bucket_tails(bucket_sizes)
        # print(last_suffix_array, "last SA")
        for char_i in np.flip(vals):
            # print(char_i)
            bucket_i = txt[char_i]
            suffix_array[bucket_tails[bucket_i]] = char_i
            bucket_tails[bucket_i] -= 1

    def build_summary_suffix_array(
        self, reduced_txt: np.ndarray, suffix_array: np.ndarray, current_names
    ):
        if current_names == reduced_txt.size:
            # print("et moi la")
            # suffix_array = np.full(reduced_txt.size + 1, 0, dtype=int)

            # suffix_array[0] = reduced_txt.size
            # for i in range(1, reduced_txt.size):
            #     suffix_array[reduced_txt[i] + 1] = i
            # #print(suffix_array, "after reduced")

            return suffix_array
        # print("jsuis la")

        return self.build_suffix_array(reduced_txt, current_names)

    def reduce(
        self,
        txt: np.ndarray,
        length: int,
        suffix_array: np.ndarray,
        suffix_types: np.ndarray,
    ):
        lms_names = np.full(length + 1, -1, dtype=int)
        current_name = 0  # size of alphabet
        counter = 1
        lms_names[suffix_array[0]] = current_name
        cur_val = suffix_array[0]
        prev_val = suffix_array[0]
        for i in range(1, length):
            if not self.is_lms_char(suffix_array[i], suffix_types):
                continue
            cur_val = suffix_array[i]
            if not self.are_lms_block_equal(
                txt, length, prev_val, cur_val, suffix_types
            ):
                current_name += 1
            prev_val = cur_val
            lms_names[cur_val] = current_name
            counter += 1

        # print(lms_names, "lms names")

        reduced_txt = np.full(counter, 0, dtype=int)
        vals = np.full(counter, 0, dtype=int)
        j = 0
        for i in range(lms_names.size):
            if lms_names[i] == -1:
                continue
            reduced_txt[j] = lms_names[i]
            vals[j] = i
            j += 1
        return (reduced_txt, vals, current_name + 1)

    def get_suffix_type(
        self,
        txt: np.ndarray,
        length: int,
    ):
        """
        Generate type S / L array
        """
        types = np.empty(length, dtype=int)
        types[length - 1] = Sais.s_type
        for i in range(length - 2, -1, -1):
            # -2 because index start at -1, and we know the last, -1 bc its not include, -1 to reverse
            if txt[i] < txt[i + 1]:
                types[i] = Sais.s_type
            elif txt[i] > txt[i + 1]:
                types[i] = Sais.l_type
            else:
                types[i] = types[i + 1]
        return types

    def get_bucket_sizes(
        self,
        txt: np.ndarray,
        length: int,
        alpha_size,
    ):
        sizes = np.full(alpha_size, 0, dtype=int)
        for i in txt:
            sizes[i] += 1

        return sizes

    def get_bucket_tails(self, bucket_size: dict):
        nb_bucket = len(bucket_size)
        tails = np.empty(nb_bucket, dtype=int)
        var = 0
        for i in range(nb_bucket):
            var += bucket_size[i]
            tails[i] = var - 1  # because not same index for bucket_size
        return tails

    def get_bucket_heads(self, bucket_size: dict):
        nb_bucket = len(bucket_size)
        head = np.empty(nb_bucket, dtype=int)
        var = 0
        for i in range(nb_bucket):
            head[i] = var  # because not same index for bucket_size
            var += bucket_size[i]
        return head

    def induction_sort_l(
        self,
        txt: np.ndarray,
        suffix_array: np.ndarray,
        suffix_types: np.ndarray,
        bucket_size: dict,
    ):
        bucket_head = self.get_bucket_heads(bucket_size)
        for i in range(suffix_array.size):
            j = suffix_array[i] - 1
            if (j < 0) or (suffix_types[j] != Sais.l_type):
                continue

            suffix_array[bucket_head[txt[j]]] = j
            bucket_head[txt[j]] += 1

    def induction_sort_s(
        self,
        txt: np.ndarray,
        suffix_array: np.ndarray,
        suffix_types: np.ndarray,
        bucket_size: dict,
    ):
        bucket_tail = self.get_bucket_tails(bucket_size)
        for i in range(suffix_array.size - 1, -1, -1):
            j = suffix_array[i] - 1
            if (j < 0) or (suffix_types[j] != Sais.s_type):
                continue
            suffix_array[bucket_tail[txt[j]]] = j
            bucket_tail[txt[j]] -= 1

    def poisition_lms_char(
        self,
        txt: np.ndarray,
        length: int,
        suffix_array: np.ndarray,
        types_list: np.ndarray,
        bucket_sizes: dict[int],
    ):
        suffix_array[0] = length - 1
        bucket_tails = self.get_bucket_tails(bucket_sizes)
        for i in range(length - 2, -1, -1):
            if not self.is_lms_char(i, types_list):
                continue
            suffix_array[bucket_tails[txt[i]]] = i
            bucket_tails[txt[i]] -= 1

    def is_lms_char(self, index: int, suffix_types: np.ndarray):
        if index == 0:
            return False
        else:
            return (
                suffix_types[index] == Sais.s_type
                and suffix_types[index - 1] == Sais.l_type
            )

    def are_lms_block_equal(
        self,
        txt: np.ndarray,
        length: int,
        prev_val,
        cur_val,
        suffix_type,
    ):
        if prev_val == length or cur_val == length:
            return False
        elif txt[prev_val] != txt[cur_val]:
            return False
        i = 1
        while i < length:
            prev_is_lms = self.is_lms_char(prev_val + i, suffix_type)
            cur_is_lms = self.is_lms_char(cur_val + i, suffix_type)
            if prev_is_lms and cur_is_lms:
                return True
            elif prev_is_lms != cur_is_lms:
                return False
            elif txt[prev_val + 1] != txt[cur_val + 1]:
                return False
            i += 1


if __name__ == "__main__":
    txt = "ACGTGCCTAGCCTACCGTGCC$"
    # txt = "GTCCCGATGTCATGTCAGGA$"
    test = Sais()
    d = tools.strToBase(txt) - 1
    a = test.build_suffix_array(d)
    # print(a)
