import getopt, sys, os


def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hi:", ["help", "input="])
    except getopt.GetoptError as err:
        print(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
    infile = None

    for opt, arg in opts:
        if opt in ("-h","--help"):
            usage()
            sys.exit()
        elif opt in ("-i", "--input"):
            if os.path.isfile(arg):
                infile = arg
        else:
            assert False, "Unhandled option"

    if infile is not None:
        splitbam(infile)


def splitbam(infile):
    acceptable_chr = ["chr" + str(x + 1) for x in range(22)]
    acceptable_chr.append("chrX")
    acceptable_chr.append("chrY")
    for chr in acceptable_chr:
        samtools_call = f"samtools view {infile} {chr} -b > {chr}_{infile}"
        print(samtools_call)
        os.system(samtools_call)

if __name__ == "__main__":
    main()
