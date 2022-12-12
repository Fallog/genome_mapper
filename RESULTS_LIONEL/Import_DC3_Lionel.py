def import_file(filename, prog):

    if prog == "DC3":
        filename = filename
    elif prog == "Map":
        filename = "Mapping_save/" + filename

    with open(filename, 'r') as f:
        # Partie DC3:
        if prog == "DC3":
            text = f.readline()
            data = text.split(" ")[:-1]
            for i in range(len(data)):
                data[i] = int(data[i])
            return (data)


if __name__ == "__main__":
    DC3_chrom1 = import_file("DC3_chrom1.txt", "DC3")
    DC3_chrom2 = import_file("DC3_chrom2.txt", "DC3")
    DC3_chrom3 = import_file("DC3_chrom3.txt", "DC3")
    DC3_chrom4 = import_file("DC3_chrom4.txt", "DC3")
    DC3_chrom5 = import_file("DC3_chrom5.txt", "DC3")
    DC3_chrom6 = import_file("DC3_chrom6.txt", "DC3")
    DC3_chrom7 = import_file("DC3_chrom7.txt", "DC3")
    DC3_chrom8 = import_file("DC3_chrom8.txt", "DC3")
    DC3_chrom9 = import_file("DC3_chrom9.txt", "DC3")
    DC3_chrom10 = import_file("DC3_chrom10.txt", "DC3")
    DC3_chrom11 = import_file("DC3_chrom11.txt", "DC3")
    DC3_chrom12 = import_file("DC3_chrom12.txt", "DC3")
    DC3_chrom13 = import_file("DC3_chrom13.txt", "DC3")
    DC3_chrom14 = import_file("DC3_chrom14.txt", "DC3")
    DC3_chrom15 = import_file("DC3_chrom15.txt", "DC3")
