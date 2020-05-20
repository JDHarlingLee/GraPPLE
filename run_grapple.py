#!/usr/bin/env python
'''Wrapper for running grapple'''

def main():
	import argparse
	
	parser = argparse.ArgumentParser(description='Arguments for main GrapPLE runner',
					 prog = 'grapple')
	
	parser.add_argument('-b', '--binary-data', required = True, type = str)
	parser.add_argument('-g', '--gene-data', required = True, type = str)
	parser.add_argument('-m', '--metadata', required = False, type = str)
	parser.add_argument('-p', '--pangenome', required = True)
	parser.add_argument('-e', '--edges', required = True)
	parser.add_argument('-t', '--threads', required = False, default = 2, type = int)
	
	parser.add_argument('-v', '--version', action = 'version', version = '%(prog)s 0.0.0')
	
	args = parser.parse_args()

	print(args)
	
if __name__ == "__main__":
	main() 

