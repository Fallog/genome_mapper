import json
import os.path


class File_handler:
    def __init__(self, file_name):
        """Create a objet with all methods to manage file with this object

        Args:
            file_name (string): basic file information
                dc3 result form : SPECIES_CHR_NUMBER
        """
        self.file_name = file_name
        self.dc3_path = "RESULTS/dc3_result"

    def export_dc3_result(self, suffix_table):
        path_file = f"{self.dc3_path}_{self.file_name}"
        if os.path.isfile(path_file):
            answer = 0
            while answer not in ("y", "n"):
                answer = input(
                    "Do you want to overwrite old dc3 result for this information ? (y/n) : "
                )
            if answer.lower() == "y":
                with open(path_file, "w") as file:
                    json.dump(suffix_table, file)
        else:
            with open(path_file, "w") as file:
                json.dump(suffix_table, file)


if __name__ == "__main__":
    a = File_handler("test_1")
    a.export_dc3_result([4, 3, 2, 6])
