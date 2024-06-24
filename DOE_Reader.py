import csv

def read_csv(path):
	configurations = []
	with open(path, newline='') as csvfile:
		reader = csv.DictReader(csvfile)
		for row in reader:
			if 'Si' in row['#consider']:
				configurations.append(row)

	return configurations