import numpy as np
import os.path


class Chromosome:
    def __init__(self, file_name, seq):
        """Create a objet with all methods to manage file with this object

        Args:
            file_name (string): basic file information
                dc3 result form : SPECIES_CHR_NUMBER
        """
        self.file_name = file_name
        self.DNA = seq
        self.DNA_dol = seq + "$"
        self.dc3_path = f"RESULTS/dc3_result_{file_name}.npy"
        self.suffix_table = None
        if os.path.isfile(self.dc3_path):
            self.import_dc3_result()

    def export_dc3_result(self, suffix_table, overwrite="y") -> None:
        """Create a .npy file contening our dc3 results

        Args:
            suffix_table (array): result of dc3 algorithm
            overwrite (str, optional):
                If "y" and file already exist, it is overwrited.
                If "n", it is not.
                If something else, it ask what the user want
                Defaults to "y".
        """
        if os.path.isfile(self.dc3_path):
            while overwrite not in ("y", "n"):
                overwrite = input(
                    "Do you want to overwrite old dc3 result for this information ? (y/n) : "
                ).lower()
            if overwrite.lower() == "y":
                self.suffix_table = suffix_table
                np.save(self.dc3_path, suffix_table)
        else:
            self.suffix_table = suffix_table
            np.save(self.dc3_path, suffix_table)

    def import_dc3_result(self, ret=True):
        """Return the suffixe table

        Returns:
            array: Suffixe table.
        """
        if os.path.isfile(self.dc3_path):
            self.suffix_table = np.load(self.dc3_path)
            return self.suffix_table
        else:
            print(f"No file at the path {self.dc3_path}.")


if __name__ == "__main__":
    a = Chromosome("test_1", "ezad")
    a.export_dc3_result([4, 3, 2, 6])
    table = a.import_dc3_result()
    print(a.suffix_table)

    print(table)
