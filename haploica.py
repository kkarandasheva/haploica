import argparse
import sys
import os
import pandas
import pybedtools

def parse_options():
    """
    Argument parser creation
    :return: options
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-i',
                        '--input-file',
                        required=True,
                        help='The file you want to be analyzed')
    parser.add_argument('-o',
                        '--output-directory',
                        required=True,
                        help='The directory into which output will be written')
    parser.add_argument('-d',
                        '--design',
                        required=True,
                        help='Amplicon panel design')
    parser.add_argument('-f',
                        '--format',
                        required=True,
                        choices=['bam', 'vcf', 'tsv'],
                        help='The input file format')
    parser.add_argument('-p',
                        '--parameter-set',
                        required=True,
                        help='Parameter dataset')
    parser_options = parser.parse_args()
    return parser_options


def err_hand_cmd(ex):
    sys.stderr.write(ex.message)
    raise SystemExit(1)


def err_hand_pipe(ex):
    raise ex


def check_file(filename, eh):
    """
    Check that file exist
    :param filename: File name to check exist
    :type filename: str
    :param eh: Error handler with single argument exception object
    :type eh: (Exception):None
    :return:
    """
    if not os.path.exists(filename):
        msg = "{0} not found\n".format(filename)
        eh(IOError(msg))


def check_path(pathname, eh):
    """
    Check that file exist
    :param pathname: Path name to check exist
    :type pathname: str
    :param eh: Error handler with single argument exception object
    :type eh: (Exception):None
    :return:
    """
    if not os.path.exists(pathname) or not os.path.isdir(pathname):
        msg = "{0} not found\n".format(pathname)
        eh(IOError(msg))


def check_extension(filename, ext, eh):
    """
    Check file extension
    :param filename:
    :param ext:
    :param eh:
    :return:
    """
    if not filename.endswith(ext):
        msg = "{0} have invalid extension\n".format(filename)
        eh(IOError(msg))


def read_parameters(filename):
    """
    Extract parameters from parameters file.
    :param filename:
    :return:
    """
    parameters = dict()
    with open(filename, 'r') as fn:
        for line in fn:
            x = line.rstrip("\r\n").split('=')
            parameters[x[0]] = eval(x[1]) if x[1].startswith('[') and x[1].endswith(']') else x[1] # TODO: FIX INJECT DANGER
    return parameters


def read_design(filename):
    df = pandas.read_csv(filename, sep = '\t', skiprows=1, names='Chromosome Start End Name Score Strand X Info'.split(' '))
    return df


def read_snv(filename, af=0.01):
    """
    Filter SNVs by AF
    :param filename:
    :param af:
    :return:
    """
    df = pandas.read_csv(filename, sep='\t')
    df = df[df.columns[[0,1,2,3,4,5,9,10,14,17]]]
    df.columns = 'Chromosome Start End Name Score Strand Ref Alt AF RS'.split(' ')
    df_filtered = df[df['AF'] > af]
    return df_filtered


def filter_snv(snv, design):
    """
    SNVs filtration by design
    :param snv:
    :param design:
    :return:
    """
    bt_snv = pybedtools.BedTool.from_dataframe(snv)
    bt_dsgn = pybedtools.BedTool.from_dataframe(design)
    result = bt_snv.intersect(bt_dsgn)
    return result


def main(infile, outdir, format, design, paramfile, errorhandler):
    """
    :param infile:
    :param outdir:
    :param format:
    :param params:
    :param errorhandler:
    :return:
    """
    check_file(infile, errorhandler)
    check_path(outdir, errorhandler)
    check_file(paramfile, errorhandler)
    check_file(infile, errorhandler)

    parameters = read_parameters(filename=paramfile)
    design = read_design(filename=design)
    snv = pandas.concat(map(read_snv, parameters['SNV'])) if isinstance(parameters['SNV'], list) else read_snv(parameters['SNV'])
    snv = filter_snv(snv=snv, design=design)

if __name__ == '__main__':
    options = parse_options()
    main(infile=options.input_file,
         outdir=options.output_directory,
         format=options.format,
         design=options.design,
         paramfile=options.parameter_set,
         errorhandler=err_hand_cmd)
    pass