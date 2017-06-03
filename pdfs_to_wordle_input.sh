## Place all of the PDFs of the papers that you want to include in a folder called
## 'pdf' within the current working directory. They will be converted to text
## inside the 'txt' directory.

mkdir -p txt
(cd pdf && batch_run.pl -0 "pdftotext #d ../txt/#d.txt")
cat txt/* > concatenated.txt

#remove some common words
perl -p -e 's/\W(et|al|figure|figures|PMID|Fig|Table|also|doi|article|thus|set)\W//gi' concatenated.txt > concatenated.clean.txt

#put plurals with singulars
perl -p -e 's/mutations/mutation/gi' concatenated.clean.txt > concatenated.clean.pluralized.txt
perl -p -e 's/populations/population/gi' concatenated.clean.pluralized.txt > concatenated.clean.pluralized.2.txt

## Now the file is ready for input into a website like
# http://www.wordle.net/
