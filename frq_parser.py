import sys

input = sys.argv[1]
output = sys.argv[2]

def read_vcf(input,output):
    with open(input, 'r') as i:
        with open(output,'w') as out:
            for line in i:
                lineOut=''
                if(line[0]=='#'):
                    line=line[2:]
                    header=line.split('\t')
                    for j in range(len(header)):
                        index = header[j].find(']')
                        colName = header[j][index + 1:]
                        if j==6:
                            lineOut=lineOut+'GENE'+'\t'+'FEATURE'+'\t'
                        elif j==len(header)-1:
                            lineOut=lineOut+colName
                        else:
                            lineOut=lineOut+colName+'\t'
                    out.write(lineOut)
                else:
                    line = line.split('\t')
                    #lineOut=''
                    for i in range(len(line)):
                        if i==0:
                            lineOut=line[i]
                        elif i != 6:
                            lineOut=lineOut+'\t'+line[i]
                        else:
                            split=line[i].split('|')
                            lineOut=lineOut+'\t'+split[0][1:]+ '\t' + split[5]
                    out.write(lineOut)

read_vcf(input, output)
