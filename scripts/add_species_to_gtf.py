__author__ = 'luca.cozzuto@crg.eu'
# -*- coding utf-8 -*-

#MODULES
import sys
import re
import optparse

#BODY FUNTIONS
def options_arg():
    usage = "usage: %prog -c <config file> -g <GTF FILE> -w <OUTPUT GTF FILE>"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option('-c', '--config-file', help='Tab separated file with species_id and species name', dest="confile" )
    parser.add_option('-g', '--gtf-file', help='GTF File (Annotation)', dest="gtffile" )
    parser.add_option('-s', '--id-size', help='id size (default 3)', dest="idisize", default=3 )
    parser.add_option('-w', '--write', help='ouput GTF File', dest="ogtf" )
    (opts,args) = parser.parse_args()
    if opts.confile and opts.gtffile and opts.ogtf:pass
    else: parser.print_help()
    return (opts)

def parse_conf(confile):
	configs = {}
	config_handle = open(confile,"r") 
	for line in config_handle.readlines():
		vals = line.split("\t")
		configs[vals[0]] = vals[1].strip()
	return (configs)

def get_attribs(string):
	values = {}
	attrs = string.strip().strip(";").split(";")
	for attr in attrs:
		rawvals = attr.strip().split(" ")	
		values[rawvals[0]] = rawvals[1]
	return (values)

def dict_to_attr_string(attrs):
	attr_string = ""
	for attrkey, attrvalue in attrs.iteritems():
		attr_string += attrkey + " " + attrvalue + "; "
	return (attr_string)
	
def fix_GTF(gtffile, idsize, species, outfile):
	outgtf = open(outfile,"w") 
	gtf_handle = open(gtffile,"r") 
	for line in gtf_handle.readlines():
		new_gtf_line = ""
		if (line[:1]!="#"):
			vals = line.split("\t")
			speciesid = vals[0][:idsize]
			species_name = species[speciesid]
			attribs = get_attribs(vals[8])
			gene_name = species_name + "_" + attribs["gene_name"].strip("\"")
			attribs["gene_name"] = "\"" + gene_name + "\""
			attr_string = dict_to_attr_string(attribs) 
			vals[8] = attr_string
			sep = "\t"
			new_gtf_line = sep.join(vals) + "\n"
		else:
			new_gtf_line = line.strip() + "\n"
		#print new_gtf_line
		outgtf.write(new_gtf_line)
	outgtf.close() 

		

     
#Calling
try:
    opts = options_arg()
    if opts.confile and opts.gtffile and opts.ogtf:
		confs = parse_conf(opts.confile)
		fix_GTF(opts.gtffile, opts.idisize, confs, opts.ogtf)
    else:
        parser.print_help()
except:pass


