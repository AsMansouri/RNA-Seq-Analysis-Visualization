import os.path
from optparse import OptionParser

def read_input(f):
    file = open(f, 'r')
    lines = list(file.read().splitlines())
    adapt = []
    k = 0
    for i in range(len(lines)):
        if len(list(lines[i])) != 0:
            if lines[i][:17] == ">>Overrepresented":
                i = i + 1
                while lines[i+1][:12] != ">>END_MODULE":
                    i = i + 1
                    k = k + 1
                    temp = lines[i].split("\t")
                    #print(">Adapt" + str(k))
                    #print(temp[0])
                    adapt.append(temp[0])
    file.close()
    return adapt

if __name__ == '__main__':

    optparser = OptionParser()
    optparser.add_option('-i', '--inputFile',
                         dest='input',
                         help='filename in .arff format',
                         default="FastQC")
    optparser.add_option('-o', '--outputFile',
                         dest='output',
                         help='The output file name',
                         default="Adapters",
                         type='str')


    (options, args) = optparser.parse_args()

    inputfile = options.input
    outputname = options.output
    file = inputfile
    Adapts = read_input(file)

    Outfile = outputname + ".fa"

    if os.path.isfile(Outfile):
        out_file = open(Outfile, 'r')
        lines = list(out_file.read().splitlines())
        k = len(lines)/2
        for a in range(len(Adapts)-1,-1,-1):
            for line in lines:
                if Adapts[a] in line:
                    del Adapts[a]
                    break
        out_file.close()
        out_file = open(Outfile, 'a')
    else:
        out_file = open(Outfile, 'w')
        k = 0

    if len(Adapts) > 0:
        for i in range(len(Adapts)):
            out_file.write('>Adapt' + str(i+k) + '\n')
            out_file.write(Adapts[i] + '\n')
    elif k != 0:
        print('No New Adapter Found ... ')
    else:
        out_file.write('>Adapt0 \n')
        out_file.write('NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n')

    out_file.close()
