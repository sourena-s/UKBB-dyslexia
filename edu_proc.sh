#for NVQ individuals with no record of age they left FTE, we use the mean of NVQ individuals who have this value: equal to 12.3 years of education
#for NVQ individuals with FTE > 25 years, we use the cap of 20 years of education (i.e. equal to university degree).
#For individuals who have no records, or chose 'prefer not to answer: -3', we designate the minimum: 7 years of education

#the input data table below includes UK Biobank fields: qualifications (data-field #6138) and age completed full time education (#845)

cat UKB-fields-data.table |awk -F"," '{q="0";if ($2~/\|-7\|/) q="none"; if ($2~/\|4\|/ || $2~/\|3\|/) q="cse"; if ($2~/\|2\|/) q="a"; if ($2~/\|6\|/) q="other"; if ($2~/\|1\|/) q="uni"; if (q=="0" && $2~/\|5\|/) q="nvq"; print $0","q;}' \
|awk -F, '{edu_years=0;if ($6=="uni") {edu_years=20;}else if ($6=="other") {edu_years=15;} else if ($6=="a") {edu_years=13;} else if ($6=="cse") {edu_years=10;}else if ($6=="none") {edu_years=7;} else if ($6=="nvq") {counter=0;edu_years=0;if ($3>4) {edu_years=$3;counter=counter+1;} if ($4>4){ edu_years=edu_years+$4;counter=counter+1;}  if ($5>4){ edu_years=edu_years+$5;counter=counter+1;};if (counter>0) {edu_years=(edu_years/counter) -5 ; if (edu_years>20) edu_years=20; } else {edu_years=12.3;}  }; if (edu_years<7) edu_years=7;    print $0","edu_years;} ' > education.txt




