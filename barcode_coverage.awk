#!/bin/awk -f
BEGIN { count = 1 }

{
	if (count == 2) { 
		prefix = substr($0, 0, 12)
		
		if (prefix in barcodes) {
			barcodes[prefix] += 1

		}
		else {
			barcodes[prefix] = 1
			
		}

	}
	else if (count == 4) {
		count = 0
	}

	count++

}

END {

	for (barcode in barcodes) {
		print (barcode, barcodes[barcode])

	}

 }
