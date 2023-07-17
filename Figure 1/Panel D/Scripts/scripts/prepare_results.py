import sys
import xlsxwriter
import pandas as pd

def add_styles(workbook):
    # formats
    red_cell = workbook.add_format({'bg_color': '#FFC7CE',
                                   'font_color': '#9C0006'})

    green_cell = workbook.add_format({'bg_color': '#C6EFCE',
                                      'font_color': '#006100'})

    yellow_cell = workbook.add_format({'bg_color': 'yellow',
                                       'bold': True})

    red_font_cell =  workbook.add_format({'font_color': 'red',
                                          'bold': True})

    cell_formats = {"red_cell":red_cell,
                    "green_cell":green_cell,
                    "yellow_cell": yellow_cell,
                    "red_font_cell":red_font_cell}

    return cell_formats



def write_mut_frequency(worksheet, mut_frequency, significance, sliding_window, styles):

    # header
    worksheet.write("A1", "Input: change cutoff if needed")
    worksheet.write("A2", "significance")
    worksheet.write("A3", significance, styles["yellow_cell"])

    worksheet.write("B2", "Codon")
    
    worksheet.write("C1", "Input: paste mutation counts")
    worksheet.write("C2", "Mutations in codon")
    
    worksheet.write("D2", "single-residue hotspot status")    
    worksheet.write("E2", "Sum of mut")
    
    worksheet.write("F1", "Sliding window calculations")
    worksheet.write("F2", "Count of mut")
    
    worksheet.write("G2", "sliding window hotspot status")    
    worksheet.write("H2", "-lg(prob)")

    # mutation frequency
    worksheet.write_column('B3', mut_frequency["Codon"])
    worksheet.write_column('C3', mut_frequency["Muts"])

    # calculate formulas
    for i in range(mut_frequency.shape[0]):
        worksheet.write_dynamic_array_formula(f"D{i+3}", f'=IFERROR(IF(C{i+3}>=$O$37,"hotspot","")," ")')
        
    for i in range(3 + sliding_window // 2, 3 + mut_frequency.shape[0] - sliding_window // 2):
        worksheet.write_dynamic_array_formula(f"E{i}", f'=SUMIF(C{i - sliding_window // 2}:C{i + sliding_window // 2},"<"&$O$37)')
        worksheet.write_dynamic_array_formula(f"F{i}", f'=COUNTIF(C{i - sliding_window // 2}:C{i + sliding_window // 2},"<"&$O$37)')
        worksheet.write_dynamic_array_formula(f"G{i}", f'=IFERROR(IF(H{i}>=(-LOG10($A$3)),"sw_htsp","")," ")')
        worksheet.write_dynamic_array_formula(f"H{i}", f'=IF(E{i}=0, 0,IFERROR(-LOG10(VLOOKUP(E{i},J:AO,5+F{i},FALSE)),-LOG10(VLOOKUP(MAX(J:J),J:AO,{sliding_window}+F{i},FALSE))))')

    for i in range(sliding_window // 2):
        worksheet.write_row(f'E{i+3}', ["=NA()"] * 4)
        worksheet.write_row(f'E{3 + mut_frequency.shape[0] - sliding_window // 2 + i}', ["=NA()"] * 4)


def write_mut_codon_frequency(worksheet):

    worksheet.write('J2', "Number of mutations in codon")
    worksheet.write('K2', "Count of codons with that many mutations")
    worksheet.write('L2', "Count of codons with less than that many mut")
    worksheet.write('M2', "Number of mut in those codons combined")

    n_muts = range(1,31)
    worksheet.write_column('J3', n_muts)
    
    for i in n_muts:
        worksheet.write(f"K{i+2}", f'=COUNTIF(C:C,J{i+2})')
        worksheet.write(f"L{i+2}", f'=COUNTIF(C:C,"<="&J{i+2})')
        worksheet.write(f"M{i+2}", f'=SUMIF(C:C,"<="&J{i+2})')

    worksheet.conditional_format(f'J3:J32',
                                 {'type': 'formula', 'criteria': '>', 'value': f'$O$37', 'format': styles["red_cell"]})


def write_single_resi_lh_threshold(worksheet):

    worksheet.write(f"K42", f'single-residue hotspots count:')
    worksheet.write(f"K43", f'=COUNTIF(C:C,">="&O37) & " codons"')
    worksheet.write(f"L43", f'="with >= "&O37&" alterations"')

    worksheet.write(f"J45", f'Count of non-hotspot codons')
    worksheet.write(f"J46", f'sum of mut in non-hotspot codons')

    worksheet.write(f"M45", f'=COUNTIF(C:C,"<"&$O$37)')
    worksheet.write(f"M46", f'=SUMIF(C:C,"<"&$O$37)')


def write_sliding_windows(worksheet, sliding_window, styles):

    n_muts = range(1,31)

    worksheet.write(f"O1", "Sliding window size, codons")
    worksheet.write(f"O2", "probability for single-residue")

    worksheet.write(f"N36", "max matching prob")
    worksheet.write(f"N37", "corresponding cutoff")
    worksheet.write(f"O35", "probabliities and cutoff values for single-residue and sliding window")

    worksheet.write_dynamic_array_formula(f"O36", f'=MAX(((O3:O32)<$A3)*(O3:O32))')
    worksheet.write_dynamic_array_formula(f"O37", f'=OFFSET($J3,MATCH(O36,O3:O32,-1)-1,0)')


    for i in n_muts:
        worksheet.write_dynamic_array_formula(f"O{i+2}", f'=IFERROR(BINOM.DIST(M{i+2}-J{i+2},M{i+2}-1,1-1/L{i+2},TRUE),1)')

    cols = ["P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z", "AA", "AB", "AC", "AD", "AE", "AF", "AG", "AH", "AI", "AJ", "AK", "AL", "AM", "AN", "AO"]

    for i in range(sliding_window-1):
        worksheet.write(f"{cols[i]}2", i + 2)
        worksheet.write_dynamic_array_formula(f"{cols[i]}36", f'=MAX((({cols[i]}3:{cols[i]}32)<$A3)*({cols[i]}3:{cols[i]}32))')
        worksheet.write_dynamic_array_formula(f"{cols[i]}37", f'=OFFSET($J3,MATCH({cols[i]}36,{cols[i]}3:{cols[i]}32,-1)-1,0)')

        for j in n_muts:
            worksheet.write_dynamic_array_formula(f"{cols[i]}{j+2}", f'=IFERROR(BINOM.DIST($M$46-$J{j+2},$M$46-1,1-{cols[i]}$2/$M$45,TRUE),1)')
    
        worksheet.conditional_format(f'{cols[i]}3:{cols[i]}32',
                                     {'type': 'cell',
                                     'criteria': 'greater than',
                                     'value':    f'${cols[i]}$36',
                                     'format': styles["red_cell"]})

        worksheet.conditional_format(f'{cols[i]}3:{cols[i]}32',
                                     {'type': 'cell',
                                     'criteria': 'less than or equal to',
                                     'value':    f'${cols[i]}$36',
                                     'format': styles["green_cell"]})

    worksheet.conditional_format(f'O3:O32',
                                 {'type': 'cell',
                                 'criteria': 'greater than',
                                 'value':    f'$O$36',
                                 'format': styles["red_cell"]})

    worksheet.conditional_format(f'O3:O32',
                                 {'type': 'cell',
                                 'criteria': 'less than or equal to',
                                 'value':    f'$O$36',
                                 'format': styles["green_cell"]})


def write_chart(worksheet, chart):
    
    chart.add_series({'values': '=Sheet1!$H$5:$H$191'})
    worksheet.insert_chart('J48', chart)


def write_3D_information(worksheet, _3d_hotspots, _3d_significance):
    data = _3d_hotspots[["Codon_1", "AA_1", "Codon_2", "AA_2", "Sum_in_group", "Count_in_group", "p_test"]]
    data = data[_3d_hotspots["p_test"] < _3d_significance]
    
    # mutation frequency
    for i, col in enumerate(data.columns):
        worksheet.write(0, i, col)
        worksheet.write_column(1, i, data[col])

    i += 1
    worksheet.write(0, i,   "Codon1_LH")
    worksheet.write(0, i+1, "Codon2_LH")
    worksheet.write(0, i+2, "Codon1_SW")
    worksheet.write(0, i+3, "Codon2_SW")


    for j in range(data.shape[0]):
        worksheet.write_dynamic_array_formula(first_row=j+1, first_col=i,   last_row=j+1, last_col=i,   formula=f'=IFNA(VLOOKUP(A{j+2},Sheet1!B:H,3,FALSE), "")')
        worksheet.write_dynamic_array_formula(first_row=j+1, first_col=i+1, last_row=j+1, last_col=i+1, formula=f'=IFNA(VLOOKUP(C{j+2},Sheet1!B:H,3,FALSE), "")')
        worksheet.write_dynamic_array_formula(first_row=j+1, first_col=i+2, last_row=j+1, last_col=i+2, formula=f'=IFNA(VLOOKUP(A{j+2},Sheet1!B:H,6,FALSE), "")')
        worksheet.write_dynamic_array_formula(first_row=j+1, first_col=i+3, last_row=j+1, last_col=i+3, formula=f'=IFNA(VLOOKUP(C{j+2},Sheet1!B:H,6,FALSE), "")')

if __name__ == '__main__':

    # intro thresholds
    _3d_significance = 0.05
    lh_significance = 0.005
    sliding_window = 5

    if sliding_window % 2 == 0:
        raise("Sliding window size should be odd")

    if sliding_window > 27:
        raise("Max size of sliding window is 27")    

    # read mut file
    mut_frequency = pd.read_csv(sys.argv[1], sep="\t")
    _3d_hotspots = pd.read_csv(sys.argv[2], sep="\t")

    # create output file
    workbook = xlsxwriter.Workbook(sys.argv[3], {'use_future_functions': True})
    styles = add_styles(workbook)

    worksheet = workbook.add_worksheet()
    write_mut_frequency(worksheet, mut_frequency, lh_significance, sliding_window, styles)
    
    write_mut_codon_frequency(worksheet)
    write_single_resi_lh_threshold(worksheet)
    chart = workbook.add_chart({'type': 'column'})
    write_sliding_windows(worksheet, sliding_window, styles)
    write_chart(worksheet, chart)

    worksheet = workbook.add_worksheet()
    write_3D_information(worksheet, _3d_hotspots, _3d_significance)


    workbook.close()



