import numpy as np
import os.path


class Chromosome:
    def __init__(self, file_name: str, seq: str, id: str):
        """Create a objet with all methods to manage file with this object

        Args:
            file_name (string): basic file information
                dc3 result form : SPECIES_CHR_NUMBER
            seq (string): DNA sequence of chromosone
            id (string): ID of the chromosome
        """
        self.id = id
        self.file_name = file_name
        self.DNA = seq
        self.DNA_dol = seq + "$"
        self.dc3_path = f"RESULTS/dc3_result_{file_name}.npy"
        self.bwt_path = f"RESULTS/bwt_result_{file_name}.txt"
        self.rank_mat_path = f"RESULTS/rank_mat_result_{file_name}.npy"
        self.suffix_table = None
        self.bwt = None
        self.rank_mat = None
        if os.path.isfile(self.dc3_path):
            self.import_dc3_result()
        if os.path.isfile(self.bwt_path):
            self.import_bwt_result()
        if os.path.isfile(self.rank_mat_path):
            self.import_rank_mat_result()

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
                if suffix_table is None:
                    np.save(self.dc3_path, self.suffix_table)
                else:
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
            if ret:
                self.suffix_table = np.load(self.dc3_path)
                return self.suffix_table
            else:
                self.suffix_table = np.load(self.dc3_path)
        else:
            print(f"No file at the path {self.dc3_path}.")

    def export_bwt_result(self, bwt=None, overwrite="y") -> None:
        """Create a .npy file contening our bwt results

        Args:
            suffix_table (array): result of bwt algorithm
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
                if bwt is None:
                    with open(self.bwt_path, "w") as f:
                        f.write(self.bwt)
                else:
                    self.bwt = bwt
                    with open(self.bwt_path, "w") as f:
                        f.write(self.bwt)
        else:
            self.bwt = bwt
            with open(self.bwt_path, "w") as f:
                f.write(self.bwt)

    def import_bwt_result(self, ret=True):
        """Return the bwt

        Returns:
            array: bwt
        """
        if os.path.isfile(self.bwt_path):
            if ret:
                with open(self.bwt_path, "r") as f:
                    self.bwt = f.read()
                return self.bwt
            else:
                with open(self.bwt_path, "r") as f:
                    self.bwt = f.read()
        else:
            print(f"No file at the path {self.bwt_path}.")

    def export_rank_matrix_result(self, rank_mat=None, overwrite="y") -> None:
        """Create a .npy file contening our bwt results

        Args:
            suffix_table (array): result of bwt algorithm
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
                if rank_mat is None:
                    np.save(self.rank_mat_path, rank_mat)

                else:
                    self.rank_mat = rank_mat
                    np.save(self.rank_mat_path, rank_mat)
        else:
            self.rank_mat = rank_mat
            np.save(self.rank_mat_path, rank_mat)

    def import_rank_mat_result(self, ret=True):
        """Return the suffixe table

        Returns:
            array: Suffixe table.
        """
        if os.path.isfile(self.rank_mat_path):
            if ret:
                self.rank_mat = np.load(self.rank_mat_path, allow_pickle=True).item()
                return self.rank_mat
            else:
                self.rank_mat = np.load(self.rank_mat_path, allow_pickle=True).item()
        else:
            print(f"No file at the path {self.rank_mat_path}.")


if __name__ == "__main__":
    a = Chromosome("test_1", "ezad", "4")
    a.export_dc3_result([4, 3, 2, 6])
    table = a.import_dc3_result()
    # a.export_rank_matrix_result({"a": 21, "b": 55})
    a.import_rank_mat_result()
    print(a.bwt, "bwt")
    print(a.suffix_table)
    print(a.rank_mat)

    print(table)
