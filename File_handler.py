import numpy as np
import os.path


class File_handler:
    def __init__(self, file_name):
        """Create a objet with all methods to manage file with this object

        Args:
            file_name (string): basic file information
                dc3 result form : SPECIES_CHR_NUMBER
        """
        self.file_name = file_name
        self.dc3_path = f"RESULTS/dc3_result_{file_name}.npy"

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
            overwrite = 0
            while overwrite not in ("y", "n"):
                overwrite = input(
                    "Do you want to overwrite old dc3 result for this information ? (y/n) : "
                ).lower()
            if overwrite.lower() == "y":
                np.save(self.dc3_path, suffix_table)
        else:
            np.save(self.dc3_path, suffix_table)

    def import_dc3_result(self):
        """Return the suffixe table

        Returns:
            array: Suffixe table.
        """
        if os.path.isfile(self.dc3_path):
            return np.load(self.dc3_path)
        else:
            print(f"No file at the path {self.dc3_path}.")


if __name__ == "__main__":
    a = File_handler("test_1")
    a.export_dc3_result([4, 3, 2, 6])
    table = a.import_dc3_result()
    print(table)
