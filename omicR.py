import tkinter as tk
from tkinter import *
from tkinter import filedialog
from tkinter import messagebox
from Bio import Entrez
import os

doubleBackSlash = r'/ '[0]


class Script1:

    def __init__(self):
        self.a_InputCSVFileName = ""
        self.b_OutputPATH = ""
        self.c_OutputFileName = "your_original_file"
        self.d_Row_Where_header_starts = 0
        self.e_rowWhereDataStarts = 0
        self.f_SEQUENCE_COLUMN = 0
        #self.g_AVG_READS = 0
        #self.i_Nbases = 0

    def Add_terminal_path_slash(self, path):
        path_with_slash = path + "\\"
        # assert isinstance(path_with_slash, object)
        return path_with_slash

    def main(self):

        # make variables
        index = 1
        index_2 = 1

        Output_path_and_name = self.b_OutputPATH + self.c_OutputFileName + ".fasta"
        Output_extra_file = self.b_OutputPATH + self.c_OutputFileName + "_InputFile_with_unique_ID.txt"
        inputFile = open(self.a_InputCSVFileName, errors='ignore')
        OutputFile = open(Output_path_and_name, 'w')
        Output_EXTRA = open(Output_extra_file, 'w')

        # MAKING FASTA FILE

        # Read through and skip header row

        for c in range(0, self.d_Row_Where_header_starts):
            HeaderTemp = inputFile.readline()

        # Reading through the rows and breaking at the end of the data

        tempstring = "temp"
        while tempstring:
            tempstring = inputFile.readline()
            if tempstring == "":
                break
            templine = tempstring.splitlines()
            x = templine[0]
            rowlist = x.split(",")
            #SeqID = rowlist[self.g_AVG_READS]
            TrimmedSequence = rowlist[self.f_SEQUENCE_COLUMN]
            #NBases = rowlist[self.i_Nbases]

            OutputRows = ">" + str(index) + '\n' + TrimmedSequence + '\n'
            #OutputRows = ">" + str(index) + "_" + SeqID + "_" + NBases + '\n' + TrimmedSequence + '\n'
            index += 1
            OutputFile.write(OutputRows)

        inputFile.close()
        OutputFile.close()

        # MAKING INPUT FILE FOR FILTERING AFTER BLAST

        # READ AND WRITE AGAIN THE INPUT FILE FOR THIS SCRIPT
        # Original File
        # Reading and writing headers

        inputFile_2 = open(self.a_InputCSVFileName, errors='ignore')
        AllHeadersJoined_inputFile_2 = ""

        for c in range(0, self.d_Row_Where_header_starts):
            headerTemp = inputFile_2.readline()

        headerLine = headerTemp.splitlines()
        y = headerLine[0]
        headerList = y.split(",")
        header_tab_delimited = ""
        for j in range(0, (len(headerList) - 1)):
            header_tab_delimited += headerList[j] + '\t'
        header_tab_delimited += headerList[(len(headerList) - 1)]
        AllHeadersJoined_inputFile_2 += ("FastaFileID" + '\t' + header_tab_delimited + '\n')  # headerList[f_SEQUENCE_COLUMN]
        Output_EXTRA.write(AllHeadersJoined_inputFile_2)

        # Original File
        # Reading through the rows and breaking at the end of the data. Writing it into a
        # new document and adding an extra column as Fasta File ID.

        tempstring = "temp"
        while tempstring:
            tempstring = inputFile_2.readline()
            if tempstring == "":
                break
            templine = tempstring.splitlines()
            x = templine[0]
            rowlist = x.split(",")
            data_tab_delimited = ""
            for i in range(0, (len(rowlist) - 1)):
                data_tab_delimited += rowlist[i] + '\t'
            data_tab_delimited += rowlist[(len(rowlist) - 1)]
            # SequenceID = rowlist[f_SEQUENCE_COLUMN]
            FastaFileID = (str(index_2))
            index_2 += 1
            data = (FastaFileID + '\t' + data_tab_delimited + '\n')
            Output_EXTRA.write(data)

        inputFile_2.close()
        Output_EXTRA.close()


class Script2:
    def __init__(self):
        self.email_DG = ""
        self.genomeAccessions_DG = ""
        self.OutputFilePath_DG = ""
        self.fileName_DG = ""

    def main(self):

        # make variables
        Entrez.email = self.email_DG

        def get_sequences_from_ID_list_line_by_line(ids):
            print(ids)

            DirectoryPath = self.OutputFilePath_DG + self.fileName_DG
            if not os.path.exists(DirectoryPath):
                os.makedirs(DirectoryPath)

            NameOfMyFile = DirectoryPath + '/' + self.fileName_DG + ".fasta"
            file_DG = open(NameOfMyFile, 'w')
            counter = 1

            for seq_id in ids:
                handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="fasta", retmode="text")
                # Read Data
                AllLines = handle.readlines()

                # PRINT AND WRITE LANE 0
                NameOfGenome_Line0 = AllLines[0].splitlines()
                print(NameOfGenome_Line0)

                str0 = ''.join(NameOfGenome_Line0)
                file_DG.write(str0)
                file_DG.write('\n')

                # Create a loop to read all rows in a file

                genome_without_header = AllLines[1:]
                listLength = len(genome_without_header)
                # print(listLength)

                complete_genome_string = ""

                for x in range(0, listLength):
                    tempList = genome_without_header[x].splitlines()
                    tempString = tempList[0]
                    complete_genome_string += tempString
                file_DG.write(complete_genome_string)
                file_DG.write('\n')

                print(counter)
                counter += 1

            file_DG.close()

        list_of_accessions = self.genomeAccessions_DG.split(',')
        get_sequences_from_ID_list_line_by_line(list_of_accessions)


class Script3:
    def __init__(self):
        self.Path_To_NCBI_BLAST_Bin_Directory = ""
        self.Path_To_Database_Fasta_File = ""
        self.Data_Base_Type = ""

    def main(self):
        # CREATE DATABASE FOR RUNNING BLAST IN WINDOWS
        CreateDataBase = self.Path_To_NCBI_BLAST_Bin_Directory + "makeblastdb -in " + self.Path_To_Database_Fasta_File + " -dbtype " + self.Data_Base_Type
        print(CreateDataBase)
        os.system(CreateDataBase)


class Script4:
    def __init__(self):
        self.x_Path_to_NCBI_Directory_BF = ""
        self.y_DC_MegaBlast_BF = 0
        self.a_Data_Base_fasta = ""
        self.b_Query_fasta_file = ""
        self.c_Output_Path_ = ""
        # BLAST PARAMETERS
        self.d_Output_file_name = "_BLAST"
        self.e_word_size = "20"
        self.f_Percentage_identity = "70"
        self.g_number_of_threads = "4"
        self.i_OutputFormat = "6"
        # FIlTERING PARAMETERS
        self.j_Percentage_overlap = "0.8"
        self.k_bitscore = "50"
        self.l_InputFile_with_unique_ID = ""

    def main(self):

        # make variables
        Task_megaBlast = ""
        if self.y_DC_MegaBlast_BF == 1:
            Task_megaBlast = " -task dc-megablast "
            print(Task_megaBlast)

        CommandLine_BF = (self.x_Path_to_NCBI_Directory_BF + "blastn " + Task_megaBlast
                       + " -db " + self.a_Data_Base_fasta + " -query "
                       + self.b_Query_fasta_file + " -out " + self.c_Output_Path_ + self.d_Output_file_name + "BLAST.txt"
                        + " -word_size " + self.e_word_size + " -perc_identity " + self.f_Percentage_identity
                         + " -num_threads " + self.g_number_of_threads + " -outfmt " + '"' + self.i_OutputFormat +
                        ' qseqid sacc stitle qseq sseq nident mismatch pident length evalue bitscore qstart qend sstart send gapopen gaps qlen slen"')

        print(CommandLine_BF)
        os.system(CommandLine_BF)

        ###################################################################################################################################
        #####################################################PART 2

        # FILTERING BLAST OUTPUTFILE

        # BLAST FILTERING PARAMETRES

        qseqid = 0
        sacc = 1
        stitle = 2
        qseq = 3
        sseq = 4
        nident = 5
        mismatch = 6
        pident = 7
        length = 8
        evalue = 9
        bitscore = 10
        qstart = 11
        qend = 12
        sstart = 13
        send = 14
        gapopen = 15
        gaps = 16
        qlen = 17
        slen = 18
        PercentageOverlapINT = 19

        BLAST_OUTPUT_FILE_BF = self.c_Output_Path_ + self.d_Output_file_name + "BLAST.txt"

        file_BF = open(BLAST_OUTPUT_FILE_BF, 'r')
        filtered_file_BF = self.c_Output_Path_ + self.d_Output_file_name + "_filtered.txt"
        filtered_files_BF = open(filtered_file_BF, "w+")

        # headers
        AllHeadersFromFilteredFile_BF = ""
        AllHeadersFromFilteredFile_BF = ("qseqid" + '\t' + "sacc" + '\t' + "stitle" + '\t' +
                                      "qseq" + '\t' + "sseq" + '\t' + "nident" + '\t' + "mismatch" + '\t' +
                                      "pident" + '\t' + "length" + '\t' + "evalue" + '\t' +
                                      "bitscore" + '\t' + "qstart" + '\t' + "qend" + '\t' + "sstart" + '\t' +
                                      "send" + '\t' + "gapopen" + '\t' + "gaps" + '\t' +
                                      "qlen" + '\t' + "slen" + '\t' + "PercentageOverlap" + '\n')
        filtered_files_BF.write(AllHeadersFromFilteredFile_BF)

        # Reading files

        tempstring = "temp"
        while tempstring:
            tempstring = file_BF.readline()
            if tempstring == "":
                break
            templine = tempstring.splitlines()

            x = templine[0]
            rowlist = x.split('\t')

            columns = (rowlist[qseqid] + '\t' + rowlist[sacc] + '\t' + rowlist[stitle] + '\t' +
                       rowlist[qseq] + '\t' + rowlist[sseq] + '\t' + rowlist[nident] + '\t' + rowlist[mismatch] + '\t' +
                       rowlist[pident] + '\t' + rowlist[length] + '\t' + rowlist[evalue] + '\t' +
                       rowlist[bitscore] + '\t' + rowlist[qstart] + '\t' + rowlist[qend] + '\t' + rowlist[
                           sstart] + '\t' +
                       rowlist[send] + '\t' + rowlist[gapopen] + '\t' + rowlist[gaps] + '\t' +
                       rowlist[qlen] + '\t' + rowlist[slen] + '\t')

            Querylength_BF = int(rowlist[qlen])
            Length_BF = int(rowlist[length])
            SubjectLength_BF = int(rowlist[slen])

            min_length_BF = min(Querylength_BF, SubjectLength_BF)

            PercentageOverlap = (Length_BF / min_length_BF)

            rowlist.append(str(PercentageOverlap))
            columns = (rowlist[qseqid] + '\t' + rowlist[sacc] + '\t' + rowlist[stitle] + '\t' +
                       rowlist[qseq] + '\t' + rowlist[sseq] + '\t' + rowlist[nident] + '\t' + rowlist[mismatch] + '\t' +
                       rowlist[pident] + '\t' + rowlist[length] + '\t' + rowlist[evalue] + '\t' +
                       rowlist[bitscore] + '\t' + rowlist[qstart] + '\t' + rowlist[qend] + '\t' + rowlist[
                           sstart] + '\t' +
                       rowlist[send] + '\t' + rowlist[gapopen] + '\t' + rowlist[gaps] + '\t' +
                       rowlist[qlen] + '\t' + rowlist[slen] + '\t' + rowlist[PercentageOverlapINT] + '\n')

            # FILTERING STEP 1 <<<<< DEFAULT "Percentage overlap >80% or 0.8" >>>>> AND <<<<< DEFAULT "BitScore >50" >>>>>

            # HANDLES

            if float(rowlist[PercentageOverlapINT]) >= float(self.j_Percentage_overlap):
                if float(rowlist[bitscore]) >= int(self.k_bitscore):
                    filtered_files_BF.write(columns)

        file_BF.close()
        filtered_files_BF.close()

        # TO BE CHECKED

        filtered_files_2_BF = open(filtered_file_BF, 'r')

        ###################################################################################################################################
        #####################################################PART 3

        # FILTERING STEP 2

        filter_part2_path_and_name_BF = self.c_Output_Path_ + self.d_Output_file_name + "_sorted.txt"
        # print(filter_part2_path_and_name)
        filtered_files_part2_BF = open(filter_part2_path_and_name_BF, "w")

        # headers
        AllHeadersFromFilteredFile = ""
        AllHeadersFromFilteredFile_BF = ("qseqid" + '\t' + "sacc" + '\t' + "stitle" + '\t' +
                                      "qseq" + '\t' + "sseq" + '\t' + "nident" + '\t' + "mismatch" + '\t' +
                                      "pident" + '\t' + "length" + '\t' + "evalue" + '\t' +
                                      "bitscore" + '\t' + "qstart" + '\t' + "qend" + '\t' + "sstart" + '\t' +
                                      "send" + '\t' + "gapopen" + '\t' + "gaps" + '\t' +
                                      "qlen" + '\t' + "slen" + '\t' + "PercentageOverlap" + '\n')
        filtered_files_part2_BF.write(AllHeadersFromFilteredFile_BF)

        # Reading files

        lst_lst = []

        counter = 0

        tempstring = "temp"
        while tempstring:
            tempstring = filtered_files_2_BF.readline()
            if tempstring == "":
                break
            if counter != 0:
                templine = tempstring.splitlines()

                x = templine[0]
                rowlist_2 = x.split('\t')

                lst_lst.append(rowlist_2)

                columns = (rowlist_2[qseqid] + '\t' + rowlist_2[sacc] + '\t' + rowlist_2[stitle] + '\t' +
                           rowlist_2[qseq] + '\t' + rowlist_2[sseq] + '\t' + rowlist_2[nident] + '\t' + rowlist_2[
                               mismatch] + '\t' +
                           rowlist_2[pident] + '\t' + rowlist_2[length] + '\t' + rowlist_2[evalue] + '\t' +
                           rowlist_2[bitscore] + '\t' + rowlist_2[qstart] + '\t' + rowlist_2[qend] + '\t' + rowlist_2[
                               sstart] + '\t' +
                           rowlist_2[send] + '\t' + rowlist_2[gapopen] + '\t' + rowlist_2[gaps] + '\t' +
                           rowlist_2[qlen] + '\t' + rowlist_2[slen] + '\t' + rowlist_2[PercentageOverlapINT] + '\n')

            counter += 1

            # READ THE NEW FILE AND ENTER THE LOOP

        # SORTING

        list.sort(lst_lst, key=lambda DataRow_0: float(DataRow_0[pident]), reverse=True)
        list.sort(lst_lst, key=lambda DataRow_2: float(DataRow_2[PercentageOverlapINT]), reverse=True)
        list.sort(lst_lst, key=lambda DataRow_1: float(DataRow_1[bitscore]), reverse=True)
        list.sort(lst_lst, key=lambda DataRow_3: DataRow_3[qseqid])

        Dictionary_lst_lst = {}

        # Reading list_list
        length = len(lst_lst)
        for i in range(length):
            temp_rowlist = lst_lst[i]

            temp_rowlist_length = len(temp_rowlist)

            if temp_rowlist_length != 20:
                print("length of tem_row_list_is:")
                print(temp_rowlist_length)
                continue

            row_string_for_output = ""

            Variable_QSeqID = temp_rowlist[qseqid]

            try:
                for j in range(temp_rowlist_length - 1):
                    temp_string = temp_rowlist[j]
                    row_string_for_output += (temp_string + "\t")

                row_string_for_output += temp_rowlist[temp_rowlist_length - 1]
                row_string_for_output += "\n"

            except IndexError:
                print("Exception thrown")
            print(row_string_for_output)

            # Tuple
            TheTuple_rowlist = (Variable_QSeqID, row_string_for_output)

            if Variable_QSeqID in Dictionary_lst_lst:
                print("key already in dictionary")
            else:
                Dictionary_lst_lst[Variable_QSeqID] = row_string_for_output
                filtered_files_part2_BF.write(row_string_for_output)

        filtered_files_part2_BF.close()
        filtered_files_2_BF.close()

        ###################################################################################################################################
        #####################################################PART 4

        ###Writting_BLAST_results_back_into_original_file

        Output_extra_file_BF = self.l_InputFile_with_unique_ID
        Output_file_only_sequences_with_hits_BF = self.c_Output_Path_ + self.d_Output_file_name + "_only_sequences_with_hits.txt"
        Output_final_BLAST_File_BF = self.c_Output_Path_ + self.d_Output_file_name + "_all_sequences_with_and_without_hits.txt"
        OutputFile_BF = open(Output_final_BLAST_File_BF, 'w')
        inputFilteredBLAST_File_BF = open(filter_part2_path_and_name_BF, 'r')

        # INDEX

        qseqid = 0
        sacc = 1
        stitle = 2
        qseq = 3
        sseq = 4
        nident = 5
        mismatch = 6
        pident = 7
        length = 8
        evalue = 9
        bitscore = 10
        qstart = 11
        qend = 12
        sstart = 13
        send = 14
        gapopen = 15
        gaps = 16
        qlen = 17
        slen = 18
        PercentageOverlapINT = 19

        # Reading files LISTS

        lst_lst = []

        counter = 0
        header_temp = ""
        Complete_output = ""

        tempstring = "temp"
        while tempstring:
            tempstring = inputFilteredBLAST_File_BF.readline()
            if counter == 0:
                Split_list = tempstring.splitlines()
                header_temp = Split_list[0]

            if tempstring == "":
                break
            if counter != 0:
                templine = tempstring.splitlines()

                x = templine[0]

                rowlist_2 = x.split('\t')

                lst_lst.append(rowlist_2)

                columns = (rowlist_2[qseqid] + '\t' + rowlist_2[sacc] + '\t' + rowlist_2[stitle] + '\t' +
                           rowlist_2[qseq] + '\t' + rowlist_2[sseq] + '\t' + rowlist_2[nident] + '\t' + rowlist_2[
                               mismatch] + '\t' +
                           rowlist_2[pident] + '\t' + rowlist_2[length] + '\t' + rowlist_2[evalue] + '\t' +
                           rowlist_2[bitscore] + '\t' + rowlist_2[qstart] + '\t' + rowlist_2[qend] + '\t' + rowlist_2[
                               sstart] + '\t' +
                           rowlist_2[send] + '\t' + rowlist_2[gapopen] + '\t' + rowlist_2[gaps] + '\t' +
                           rowlist_2[qlen] + '\t' + rowlist_2[slen] + '\t' + rowlist_2[PercentageOverlapINT] + '\n')

            counter += 1

        Dictionary_lst_lst = {}

        # Reading list_list

        length = len(lst_lst)
        for i in range(length):
            temp_rowlist = lst_lst[i]

            temp_rowlist_length = len(temp_rowlist)
            if temp_rowlist_length != 20:
                continue

            row_string_for_output = ""

            Variable_QSeqID = temp_rowlist[qseqid]

            try:
                for j in range(temp_rowlist_length):
                    temp_string = temp_rowlist[j]
                    row_string_for_output += (temp_string + "\t")
                row_string_for_output += "\n"
            except IndexError:

                print("Exception thrown")

            # Tuple
            TheTuple_rowlist = (Variable_QSeqID, row_string_for_output)

            if Variable_QSeqID in Dictionary_lst_lst:
                print("key already in dictionary")
            else:
                Dictionary_lst_lst[Variable_QSeqID] = row_string_for_output
                print(row_string_for_output)

        # OPEN THE ORIGINAL MODIFIED FILE

        Original_Modified_file_BF = open(Output_extra_file_BF, 'r')
        Only_sequences_with_hits_file_BF = open(Output_file_only_sequences_with_hits_BF, 'w')

        counter2 = 0
        Header_Temp_2 = ""

        tempstring = "temp"
        while tempstring:
            tempstring = Original_Modified_file_BF.readline()

            if counter2 == 0:
                Split_list_2 = tempstring.splitlines()
                Header_Temp_2 = Split_list_2[0]
                OutputFile_BF.write(Header_Temp_2 + "\t" + header_temp + "\n")
                Only_sequences_with_hits_file_BF.write(Header_Temp_2 + "\t" + header_temp + "\n")

            if tempstring == "":
                break

            if counter2 != 0:
                templine = tempstring.splitlines()
                x = templine[0]

                rowlist = x.split('\t')
                Temp_QSeqID = rowlist[0]

                if Temp_QSeqID in Dictionary_lst_lst:
                    Corresponding_row = Dictionary_lst_lst.get(Temp_QSeqID)

                    OutputFile_BF.write(x + "\t")
                    OutputFile_BF.write(Corresponding_row)
                    Only_sequences_with_hits_file_BF.write(x + "\t")
                    Only_sequences_with_hits_file_BF.write(Corresponding_row)

                else:
                    OutputFile_BF.write(
                        x + "\t" + "\t" + "\t" + "\t" + "\t" + "\t" + "\t" + "\t" + "\t" + "\t" + "\t" + "\t" + "\t" + "\t" + "\t" + "\t" + "\t" + "\t" + "\t" + "\t" + "\n")
                    # print("not found in dictionary")

            counter2 += 1

        OutputFile_BF.write(Complete_output)

        Original_Modified_file_BF.close()
        OutputFile_BF.close()
        inputFilteredBLAST_File_BF.close()
        Only_sequences_with_hits_file_BF.close()
        Only_sequences_with_hits_file_BF.close()


class Script55:
    def __init__(self):

        self.a_BLAST_input_path_and_file_ = ""
        self.b_Output_Path_ = ""
        self.c_Output_file_name = "_BLAST"

        # FIlTERING PARAMETERS
        self.d_Percentage_overlap = "0.8"
        self.e_bitscore = "50"

    def main(self):

        # BLAST FILTERING PARAMETERS

        qseqid = 0
        sacc = 1
        stitle = 2
        qseq = 3
        sseq = 4
        nident = 5
        mismatch = 6
        pident = 7
        length = 8
        evalue = 9
        bitscore = 10
        qstart = 11
        qend = 12
        sstart = 13
        send = 14
        gapopen = 15
        gaps = 16
        qlen = 17
        slen = 18
        PercentageOverlapINT = 19

        BLAST_OUTPUT_FILE_F = self.a_BLAST_input_path_and_file_

        file_F = open(BLAST_OUTPUT_FILE_F, 'r')
        filtered_file_F = self.b_Output_Path_ + self.c_Output_file_name + "_filtered.txt"
        filtered_files_F = open(filtered_file_F, "w+")

        # headers
        AllHeadersFromFilteredFile_F = ""
        AllHeadersFromFilteredFile_F = ("qseqid" + '\t' + "sacc" + '\t' + "stitle" + '\t' +
                                      "qseq" + '\t' + "sseq" + '\t' + "nident" + '\t' + "mismatch" + '\t' +
                                      "pident" + '\t' + "length" + '\t' + "evalue" + '\t' +
                                      "bitscore" + '\t' + "qstart" + '\t' + "qend" + '\t' + "sstart" + '\t' +
                                      "send" + '\t' + "gapopen" + '\t' + "gaps" + '\t' +
                                      "qlen" + '\t' + "slen" + '\t' + "PercentageOverlap" + '\n')
        filtered_files_F .write(AllHeadersFromFilteredFile_F)

        # Reading files

        tempstring = "temp"
        while tempstring:
            tempstring = file_F.readline()
            if tempstring == "":
                break
            templine = tempstring.splitlines()

            x = templine[0]
            rowlist = x.split('\t')

            columns = (rowlist[qseqid] + '\t' + rowlist[sacc] + '\t' + rowlist[stitle] + '\t' +
                       rowlist[qseq] + '\t' + rowlist[sseq] + '\t' + rowlist[nident] + '\t' + rowlist[mismatch] + '\t' +
                       rowlist[pident] + '\t' + rowlist[length] + '\t' + rowlist[evalue] + '\t' +
                       rowlist[bitscore] + '\t' + rowlist[qstart] + '\t' + rowlist[qend] + '\t' + rowlist[
                           sstart] + '\t' +
                       rowlist[send] + '\t' + rowlist[gapopen] + '\t' + rowlist[gaps] + '\t' +
                       rowlist[qlen] + '\t' + rowlist[slen] + '\t')

            Querylength_F = int(rowlist[qlen])
            Length_F = int(rowlist[length])
            SubjectLength_F = int(rowlist[slen])

            min_length_F = min(Querylength_F, SubjectLength_F)

            PercentageOverlap_F = (Length_F / min_length_F)

            rowlist.append(str(PercentageOverlap_F))
            columns = (rowlist[qseqid] + '\t' + rowlist[sacc] + '\t' + rowlist[stitle] + '\t' +
                       rowlist[qseq] + '\t' + rowlist[sseq] + '\t' + rowlist[nident] + '\t' + rowlist[mismatch] + '\t' +
                       rowlist[pident] + '\t' + rowlist[length] + '\t' + rowlist[evalue] + '\t' +
                       rowlist[bitscore] + '\t' + rowlist[qstart] + '\t' + rowlist[qend] + '\t' + rowlist[
                           sstart] + '\t' +
                       rowlist[send] + '\t' + rowlist[gapopen] + '\t' + rowlist[gaps] + '\t' +
                       rowlist[qlen] + '\t' + rowlist[slen] + '\t' + rowlist[PercentageOverlapINT] + '\n')

            # FILTERING STEP 1 <<<<< DEFAULT "Percentage overlap >80% or 0.8" >>>>> AND <<<<< DEFAULT "BitScore >50" >>>

            # HANDLES

            if float(rowlist[PercentageOverlapINT]) >= float(self.d_Percentage_overlap):
                if float(rowlist[bitscore]) >= int(self.e_bitscore):
                    filtered_files_F.write(columns)

        file_F.close()
        filtered_files_F.close()

        # TO BE CHECKED

        filtered_files_2_F = open(filtered_file_F, 'r')

        ###################################################################################################################################
        #####################################################PART 3

        # FILTERING STEP 2

        filter_part2_path_and_name_F = self.b_Output_Path_ + self.c_Output_file_name + "_sorted.txt"
        filtered_files_part2_F = open(filter_part2_path_and_name_F, "w")

        # headers
        AllHeadersFromFilteredFile_F = ""
        AllHeadersFromFilteredFile_F = ("qseqid" + '\t' + "sacc" + '\t' + "stitle" + '\t' +
                                      "qseq" + '\t' + "sseq" + '\t' + "nident" + '\t' + "mismatch" + '\t' +
                                      "pident" + '\t' + "length" + '\t' + "evalue" + '\t' +
                                      "bitscore" + '\t' + "qstart" + '\t' + "qend" + '\t' + "sstart" + '\t' +
                                      "send" + '\t' + "gapopen" + '\t' + "gaps" + '\t' +
                                      "qlen" + '\t' + "slen" + '\t' + "PercentageOverlap" + '\n')
        filtered_files_part2_F.write(AllHeadersFromFilteredFile_F)

        # Reading files

        lst_lst = []

        counter = 0

        tempstring = "temp"
        while tempstring:
            tempstring = filtered_files_2_F.readline()
            if tempstring == "":
                break
            if counter != 0:
                templine = tempstring.splitlines()

                x = templine[0]
                rowlist_2 = x.split('\t')

                lst_lst.append(rowlist_2)

                columns = (rowlist_2[qseqid] + '\t' + rowlist_2[sacc] + '\t' + rowlist_2[stitle] + '\t' +
                           rowlist_2[qseq] + '\t' + rowlist_2[sseq] + '\t' + rowlist_2[nident] + '\t' + rowlist_2[
                               mismatch] + '\t' +
                           rowlist_2[pident] + '\t' + rowlist_2[length] + '\t' + rowlist_2[evalue] + '\t' +
                           rowlist_2[bitscore] + '\t' + rowlist_2[qstart] + '\t' + rowlist_2[qend] + '\t' + rowlist_2[
                               sstart] + '\t' +
                           rowlist_2[send] + '\t' + rowlist_2[gapopen] + '\t' + rowlist_2[gaps] + '\t' +
                           rowlist_2[qlen] + '\t' + rowlist_2[slen] + '\t' + rowlist_2[PercentageOverlapINT] + '\n')

            counter += 1

            # READ THE NEW FILE AND ENTER THE LOOP

        # SORTING

        list.sort(lst_lst, key=lambda DataRow_0: float(DataRow_0[pident]), reverse=True)
        list.sort(lst_lst, key=lambda DataRow_2: float(DataRow_2[PercentageOverlapINT]), reverse=True)
        list.sort(lst_lst, key=lambda DataRow_1: float(DataRow_1[bitscore]), reverse=True)
        list.sort(lst_lst, key=lambda DataRow_3: DataRow_3[qseqid])

        Dictionary_lst_lst = {}

        # Reading list_list
        length = len(lst_lst)
        for i in range(length):
            temp_rowlist = lst_lst[i]

            temp_rowlist_length = len(temp_rowlist)

            if temp_rowlist_length != 20:
                print("length of tem_row_list_is:")
                print(temp_rowlist_length)
                continue

            row_string_for_output = ""

            Variable_QSeqID = temp_rowlist[qseqid]

            try:
                for j in range(temp_rowlist_length - 1):
                    temp_string = temp_rowlist[j]
                    row_string_for_output += (temp_string + "\t")

                row_string_for_output += temp_rowlist[temp_rowlist_length - 1]
                row_string_for_output += "\n"

            except IndexError:
                print("Exception thrown")
            print(row_string_for_output)

            # Tuple
            TheTuple_rowlist = (Variable_QSeqID, row_string_for_output)

            if Variable_QSeqID in Dictionary_lst_lst:
                print("key already in dictionary")
            else:
                Dictionary_lst_lst[Variable_QSeqID] = row_string_for_output
                filtered_files_part2_F.write(row_string_for_output)

        filtered_files_part2_F.close()
        filtered_files_2_F.close()


class Win1(Script1, Script2, Script3, Script4, Script55):
    def __init__(self, window):
        # Initializations
        self.wind = window
        self.wind.title("omicR")
        self.wind.wm_iconbitmap('Currito.ico')
        self.wind.resizable(False, False)

        # Creating A Frame Container
        frame = LabelFrame(self.wind, text="Select what would you like to do:")
        frame.grid(row=0, column=0, columnspan=3, padx=40, pady=40)

        # Buttons
        tk.Button(frame, text="Create FASTA files and input files for BLAST / filtering",
                  command=self.new_window2).grid(row=3, columnspan=2, padx=5, pady=5, sticky=W + E)

        tk.Button(frame, text="Download Genomes",
                  command=self.new_window3).grid(row=4, columnspan=2, padx=5, pady=5, sticky=W + E)

        tk.Button(frame, text="Create Genome Database",
                  command=self.new_window4).grid(row=5, columnspan=2, padx=5, pady=5, sticky=W + E)

        tk.Button(frame, text="BLAST / filtering", command=self.new_window5
                  ).grid(row=6, columnspan=2, padx=5, pady=5, sticky=W + E)

        tk.Button(frame, text="Filtering", command=self.new_window55
                  ).grid(row=7, columnspan=2, padx=5, pady=5, sticky=W + E)

        # Instructions
        tk.Button(frame, text="Instructions", command=self.new_window6).grid(row=8, columnspan=2, padx=5, pady=5,
                                                                             sticky=W + E)
        # Close Button
        tk.Button(frame, text="Close", command=self.close_window).grid(row=10, column=1, columnspan=2, padx=5, pady=5,
                                                                       sticky=E)

    # FASTA FILES
    def new_window2(self):
        self.new_window = tk.Toplevel(self.wind)
        self.app = Win2(self.new_window)

    # Download Genomes
    def new_window3(self):
        self.new_window = tk.Toplevel(self.wind)
        self.app = Win3(self.new_window)

    # Create Genome Database
    def new_window4(self):
        self.new_window = tk.Toplevel(self.wind)
        self.app = Win4(self.new_window)

    # BLAST and Filtering
    def new_window5(self):
        self.new_window = tk.Toplevel(self.wind)
        self.app = Win5(self.new_window)

    def new_window55(self):
        self.new_window = tk.Toplevel(self.wind)
        self.app = Win55(self.new_window)

    # HELP
    def new_window6(self):
        self.new_window = tk.Toplevel(self.wind)
        self.app = Win6(self.new_window)

    def close_window(self):
        self.wind.destroy()


# FASTA FILES
class Win2(Win1):

    def __init__(self, window):
        # Initializations
        self.window = window
        self.wind = window
        self.wind.title("Create FASTA files and input files for BLAST filtering")
        self.wind.wm_iconbitmap('Currito.ico')

        # Creating a Frame Container
        frame = LabelFrame(self.wind, text="Complete the following parameters ")
        frame.grid(row=0, column=0, columnspan=3, padx=20, pady=20)

        # CSV Input
        # ROW1
        # INPUT FILE PATH
        Label(frame, text="Input file path (CSV file required): ").grid(row=1, column=0, sticky=W)
        Label(frame, text="Example: " + "C:" + doubleBackSlash +
                          doubleBackSlash + "Users/MyDocuments/Bassiana.csv ").grid(row=1, column=6, sticky=W)
        self.CSVInput = tk.Entry(frame)
        self.CSVInput.focus()
        self.CSVInput.grid(row=1, column=1, columnspan=2, ipadx=200, padx=5, pady=5, )
        # Button Select
        tk.Button(frame, text="Select", command=self.select_directory_CSV_input_path).grid(row=1, column=5,
                                                                                           sticky=W + E, padx=2, pady=2)

        # Output File Path
        # ROW 2
        Label(frame, text="Output file path: ").grid(row=2, column=0, sticky=W)
        Label(frame, text="Example: " + "C:" + doubleBackSlash +
                          doubleBackSlash + "Users/MyDocuments/ ").grid(row=2, column=6, sticky=W)
        self.OutputFilePath = tk.Entry(frame)
        self.OutputFilePath.grid(row=2, column=1, columnspan=2, ipadx=200, padx=5, pady=5, )
        # Button Select
        tk.Button(frame, text="Select", command=self.select_directory_output_path).grid(row=2, column=5, padx=2, pady=2)

        # OUTPUT FILE NAME
        # ROW 3
        Label(frame, text="Output file name: ").grid(row=3, column=0, sticky=W)
        Label(frame, text="Example: Bassiana_BLAST_Results ").grid(row=3, column=6, sticky=W)
        self.CSVOutputFileName = tk.Entry(frame)
        self.CSVOutputFileName.grid(row=3, column=1, columnspan=2, ipadx=200, padx=5, pady=5, )

        # Row where header starts
        # Row 4
        Label(frame, text="Row where header starts: ").grid(row=4, column=0, sticky=W)
        Label(frame, text="Example: 1 (Start counting rows from one *Number*").grid(row=4, column=2, sticky=W)
        self.RowWhereHeaderStarts = tk.Entry(frame)
        self.RowWhereHeaderStarts.grid(row=4, column=1, columnspan=1, ipadx=50, padx=5, pady=5, sticky=W)

        # Row where data starts
        # Row 5
        Label(frame, text="Row where data starts: ").grid(row=5, column=0, sticky=W)
        Label(frame, text="Example: 2 (Start counting rows from one *Number*)").grid(row=5, column=2, sticky=W)
        self.RowWhereDataStarts = tk.Entry(frame)
        self.RowWhereDataStarts.grid(row=5, column=1, columnspan=1, ipadx=50, padx=5, pady=5, sticky=W)

        # Column of sequences
        # Row 6
        Label(frame, text="Column of Sequences: ").grid(row=6, column=0, sticky=W)
        Label(frame, text="Example: 0 (Start counting columns from zero *Index*)").grid(row=6, column=2, sticky=W)
        self.ColumnOfSequences = tk.Entry(frame)
        self.ColumnOfSequences.grid(row=6, column=1, columnspan=1, ipadx=50, padx=5, pady=5, sticky=W)

        # Column of AVGreads
        # Row 7
        #Label(frame, text="Column of comments [Average Reads]: ").grid(row=7, column=0, sticky=W)
        #Label(frame, text="Example: 1 (Start counting columns from zero *Index*)").grid(row=7, column=2, sticky=W)
        #self.ColumnOfAVGreads = tk.Entry(frame)
        #self.ColumnOfAVGreads.grid(row=7, column=1, columnspan=1, ipadx=50, padx=5, pady=5, sticky=W)

        # Column of N Bases
        # Row 8
        #Label(frame, text="Column of comments [NBases]: ").grid(row=8, column=0, sticky=W)
        #Label(frame, text="Example: 2 (Start counting columns from zero *Index*)").grid(row=8, column=2, sticky=W)
        #Label(frame, text=" ").grid(row=9, column=6, sticky=W)
        #self.ColumnOfNbases = tk.Entry(frame)
        #self.ColumnOfNbases.grid(row=8, column=1, columnspan=1, ipadx=50, padx=5, pady=5, sticky=W)

        # WINDOW FASTA FILES

        # Button clear
        tk.Button(frame, text="Clear all", command=lambda: [self.ClearAll_Fasta_files()]).grid(row=17, column=6,
                                                                                               columnspan=1, padx=5,
                                                                                               pady=5,
                                                                                               sticky=W + E)

        # Button Run
        tk.Button(frame, text="Run", command=lambda: [self.Run_Button_FASTA_FILES()]).grid(row=18, column=6,
                                                                                           columnspan=1, padx=5, pady=5,
                                                                                           sticky=W + E)

        # BUTTON Close
        tk.Button(frame, text="Close", command=self.close_window).grid(row=19, column=6, columnspan=1, padx=5, pady=5,
                                                                       sticky=W + E)

    def Run_Button_FASTA_FILES(self):
        TextFromCSVInput = self.CSVInput.get()
        TextFromOutputPath = self.OutputFilePath.get()
        TextFromOutputFileName = self.CSVOutputFileName.get()
        TextFromRowWhereHeaderStarts = self.RowWhereHeaderStarts.get()
        TextFromRowWhereDataStarts = self.RowWhereDataStarts.get()
        TextFromSequenceColumn = self.ColumnOfSequences.get()
        #TextFromAVG_reads = self.ColumnOfAVGreads.get()
        #TextFrom_Nbases = self.ColumnOfNbases.get()

        if (len(TextFromCSVInput) != 0
                and len(TextFromOutputPath) != 0
                and len(TextFromRowWhereHeaderStarts) != 0
                and len(TextFromRowWhereDataStarts) != 0
                and len(TextFromSequenceColumn)):
            Script1.a_InputCSVFileName = TextFromCSVInput
            TextFromOutputMod = TextFromOutputPath + doubleBackSlash
            Script1.b_OutputPATH = TextFromOutputMod

            if len(TextFromOutputFileName) != 0:
                Script1.c_OutputFileName = TextFromOutputFileName

            Script1.c_OutputFileName = TextFromOutputFileName
            Script1.d_Row_Where_header_starts = int(TextFromRowWhereHeaderStarts)
            Script1.e_rowWhereDataStarts = int(TextFromRowWhereDataStarts)
            Script1.f_SEQUENCE_COLUMN = int(TextFromSequenceColumn)
            #Script1.g_AVG_READS = int(TextFromAVG_reads)
            #Script1.i_Nbases = int(TextFrom_Nbases)

            print(Script1.a_InputCSVFileName)
            print(Script1.b_OutputPATH)
            print(Script1.c_OutputFileName)
            print(Script1.d_Row_Where_header_starts)
            print(Script1.e_rowWhereDataStarts)
            print(Script1.f_SEQUENCE_COLUMN)
            #print(Script1.g_AVG_READS)
            #print(Script1.i_Nbases)

            # Output Messages

            self.message = Label(text="Running", fg='red').grid(row=19, column=2, columnspan=2, sticky=W + E, padx=5,
                                                                pady=5)

            Script1.main(self)

            self.message1 = Label(text="Completed!", fg='red').grid(row=19, column=2, columnspan=2, sticky=W + E,
                                                                    padx=5,
                                                                    pady=5)
            messagebox.showinfo('Information', "Completed!")

        else:
            messagebox.showerror("Error", "All required parameters must be filled")
        self.wind.lift()

    def select_directory_CSV_input_path(self):
        folder_selected = filedialog.askopenfilename(initialdir='/', title="Select file",
                                                     filetypes=(("CSV files", "*.csv"), ("all files", "*.*")))
        print(folder_selected)
        self.CSVInput.delete(0, END)
        self.CSVInput.insert(0, folder_selected)
        self.wind.lift()
        return

    def select_directory_output_path(self):
        Output_file_path_fasta_files = filedialog.askdirectory(initialdir='.')
        print(Output_file_path_fasta_files)
        self.OutputFilePath.delete(0, END)
        self.OutputFilePath.insert(0, Output_file_path_fasta_files)
        self.wind.lift()
        return

    def ClearAll_Fasta_files(self):
        self.CSVInput.delete(0, END)
        self.OutputFilePath.delete(0, END)
        self.CSVOutputFileName.delete(0, END)
        self.RowWhereHeaderStarts.delete(0, END)
        self.RowWhereDataStarts.delete(0, END)
        self.ColumnOfSequences.delete(0, END)
        #self.ColumnOfAVGreads.delete(0, END)
        #self.ColumnOfNbases.delete(0, END)

    def close_window(self):
        self.wind.destroy()


# Download Genomes
class Win3(Win1):
    def __init__(self, window):
        # Initializations
        self.window = window
        self.wind = window
        self.wind.title("Download genome entries from NCBI")
        self.wind.wm_iconbitmap('Currito.ico')

        # Creating a Frame Container
        frame = LabelFrame(self.wind, text="Complete the following parameters ")
        frame.grid(row=0, column=0, columnspan=3, padx=20, pady=20)

        # Write your email
        # Row1
        Label(frame, text="E-mail (required to access NCBI): ").grid(row=1, column=0, sticky=W)
        Label(frame, text="Example: DungogCitizen@nothing.com ").grid(row=1, column=6, sticky=W)
        self.WriteYourEmail = tk.Entry(frame)
        self.WriteYourEmail.focus()
        self.WriteYourEmail.grid(row=1, column=1, columnspan=2, ipadx=200, padx=5, pady=5, )

        # Write accession numbers
        # Row 2
        Label(frame, text="Accession numbers: ").grid(row=2, column=0, sticky=W)
        Label(frame, text="Example: NC_009328.1, NC_009329.1 ").grid(row=2, column=6, sticky=W)
        self.WriteAccessionNumbers = tk.Entry(frame)
        self.WriteAccessionNumbers.grid(row=2, column=1, columnspan=2, ipadx=200, padx=5, pady=5, )

        # Your output path for downloading genomes
        # Row 3
        Label(frame, text="Output path directory: ").grid(row=3, column=0, sticky=W)
        Label(frame, text="Example: " + "C:" + doubleBackSlash +
                          doubleBackSlash + "Users/MyDocuments/ ").grid(row=3, column=6, sticky=W)
        self.OutputPathForDownloadingGenomes = tk.Entry(frame)
        self.OutputPathForDownloadingGenomes.grid(row=3, column=1, columnspan=2, ipadx=200, padx=5, pady=5, )
        # Button Select
        tk.Button(frame, text="Select", command=self.Download_genomes_select_directory_output_path).grid(row=3,
                                                                                                         column=5,
                                                                                                         padx=2, pady=2)

        # Downloading genomes output file name
        # Row 4
        Label(frame, text="Output file name: ").grid(row=4, column=0, sticky=W)
        Label(frame, text="Example: Geobacillus_sp_Genome").grid(row=4, column=6, sticky=W)
        Label(frame, text="").grid(row=5, column=6, sticky=W)
        self.DownloadingGenomes_OutputFileName = tk.Entry(frame)
        self.DownloadingGenomes_OutputFileName.grid(row=4, column=1, columnspan=2, ipadx=200, padx=5, pady=5, )

        # Button clear
        tk.Button(frame, text="Clear all", command=lambda: [self.ClearAll_DownloadGENOMES()]).grid(row=17, column=6,
                                                                                                   columnspan=1, padx=5,
                                                                                                   pady=5,
                                                                                                   sticky=W + E)

        # Button Run
        tk.Button(frame, text="Run", command=self.Run_Button_Downloading_genomes).grid(row=18, column=6,
                                                                                       columnspan=1, padx=5,
                                                                                       pady=5, sticky=W + E)

        # Button Close
        tk.Button(frame, text="Close", command=self.close_window).grid(row=19, column=6, columnspan=1, padx=5, pady=5,
                                                                       sticky=W + E)

    def Run_Button_Downloading_genomes(self):
        TextFromWriteYourEmail_DownloadingGenomes = self.WriteYourEmail.get()
        TextFromWriteAccessionNumbers_DownloadingGenomes = self.WriteAccessionNumbers.get()
        TextFromOutputPathFor_DownloadingGenomes = self.OutputPathForDownloadingGenomes.get()
        TextFromOutputFileName_DownloadingGenomes = self.DownloadingGenomes_OutputFileName.get()

        if (len(TextFromWriteYourEmail_DownloadingGenomes) != 0
                and len(TextFromWriteAccessionNumbers_DownloadingGenomes) != 0
                and len(TextFromOutputPathFor_DownloadingGenomes) != 0
                and len(TextFromOutputFileName_DownloadingGenomes) != 0):

            Script2.email_DG = TextFromWriteYourEmail_DownloadingGenomes
            Script2.genomeAccessions_DG = TextFromWriteAccessionNumbers_DownloadingGenomes

            TextFromOutputMod_DownloadingGenomes = TextFromOutputPathFor_DownloadingGenomes + doubleBackSlash
            Script2.OutputFilePath_DG = TextFromOutputMod_DownloadingGenomes
            Script2.fileName_DG = TextFromOutputFileName_DownloadingGenomes

            print(Script2.email_DG)
            print(Script2.genomeAccessions_DG)
            print(Script2.OutputFilePath_DG)
            print(Script2.fileName_DG)

            # Output Messages

            self.message = Label(text="Running", fg='red').grid(row=19, column=2, columnspan=2, sticky=W + E, padx=5,
                                                                pady=5)
            Script2.main(self)
            self.message1 = Label(text="Completed!", fg='red').grid(row=19, column=2, columnspan=2, sticky=W + E,
                                                                    padx=5,
                                                                    pady=5)
            messagebox.showinfo('Information', "Completed!")
        else:
            messagebox.showerror("Error", "All required parameters must be filled")
        self.wind.lift()

    def Download_genomes_select_directory_output_path(self):
        Output_file_path_Downloading_genomes = filedialog.askdirectory(initialdir='.')
        print(Output_file_path_Downloading_genomes)
        self.OutputPathForDownloadingGenomes.delete(0, END)
        self.OutputPathForDownloadingGenomes.insert(0, Output_file_path_Downloading_genomes)
        self.wind.lift()
        return

    def ClearAll_DownloadGENOMES(self):
        self.WriteYourEmail.delete(0, END)
        self.WriteAccessionNumbers.delete(0, END)
        self.OutputPathForDownloadingGenomes.delete(0, END)
        self.DownloadingGenomes_OutputFileName.delete(0, END)

    def close_window(self):
        self.wind.destroy()


# Create Genome Database
class Win4(Win1):

    def __init__(self, window):
        # Initializations
        self.window = window
        self.wind = window
        self.wind.title("Create NCBI Database for BLASTn")
        self.wind.wm_iconbitmap(
            'Currito.ico')

        # Creating a Frame Container
        frame = LabelFrame(self.wind, text="Complete the following parameters")
        frame.grid(row=0, column=0, columnspan=3, padx=20, pady=20)

        # Select your NCBI Path
        # Row 1
        Label(frame, text="Select path to NCBI/bin directory : ").grid(row=1, column=0, sticky=W)
        Label(frame, text="Example: " + "C:" + doubleBackSlash +
                          doubleBackSlash + "NCBI" + doubleBackSlash + "blast-2.8.0+" + doubleBackSlash + "bin   ").grid(
            row=1, column=6, sticky=W)
        self.NCBIPath_to_BIN = tk.Entry(frame)
        self.NCBIPath_to_BIN.focus()
        self.NCBIPath_to_BIN.grid(row=1, column=1, columnspan=2, ipadx=200, padx=5, pady=5, )
        # Button Select
        tk.Button(frame, text="Select", command=self.Select_button_for_Create_Database_for_NCBI_BLAST).grid(row=1,
                                                                                                            column=5,
                                                                                                            padx=2,
                                                                                                            pady=2)

        # Path to Fasta file to Build DB
        # Row 2
        Label(frame, text="Select path to genome FASTA file: ").grid(row=2, column=0, sticky=W)
        Label(frame, text="Example: " + "C:" + doubleBackSlash + doubleBackSlash
                          + "Users/MyDocuments/Genome/Genome.fasta ").grid(row=2, column=6, sticky=W)
        self.Path_to_FASTA_file = tk.Entry(frame)
        self.Path_to_FASTA_file.grid(row=2, column=1, columnspan=2, ipadx=200, padx=5, pady=5, )
        # Button Select
        tk.Button(frame, text="Select", command=self.select_path_for_FASTA_file).grid(row=2, column=5, padx=2, pady=2)

        # Drop down menu
        # Database type
        # Row 3
        Label(frame, text="Click to select type of database: ").grid(row=3, column=0, sticky=W)
        Label(frame, text="Example: 'nucl' for nucleotide or 'prot' for protein").grid(row=3, column=2, sticky=W)
        Label(frame, text="").grid(row=4, column=6, sticky=W)

        Options = ["", "nucl", "prot"]
        self.clicked = StringVar()
        self.clicked.set(Options[0])

        self.DropMenu_DB_Type = OptionMenu(frame, self.clicked, *Options,)
        self.DropMenu_DB_Type.grid(row=3, column=1, columnspan=1, padx=5, pady=5, sticky=W)

        # Button clear
        tk.Button(frame, text="Clear all", command=lambda: [self.ClearAll_Create_GENOME_DB()]).grid(row=17, column=6,
                                                                                                    columnspan=1,
                                                                                                    padx=5,
                                                                                                    pady=5,
                                                                                                    sticky=W + E)
        # Button Run
        tk.Button(frame, text="Run", command=self.Run_Button_Create_Database).grid(row=18, column=6, columnspan=1,
                                                                                   padx=5, pady=5, sticky=W + E)

        # BUTTON Close
        tk.Button(frame, text="Close", command=self.close_window).grid(row=19, column=6, columnspan=1, padx=5, pady=5,
                                                                       sticky=W + E)

    def Run_Button_Create_Database(self):
        TextFromNCBI_path_to_bin_directory = self.NCBIPath_to_BIN.get()
        TextFrom_Path_to_Genome_fasta_file = self.Path_to_FASTA_file.get()
        TextFromDropDownMenu_DataBaseType = self.clicked.get()

        if (len(TextFromNCBI_path_to_bin_directory) != 0
                and len(TextFrom_Path_to_Genome_fasta_file) != 0
                and len(TextFromDropDownMenu_DataBaseType) !=0):


            TextFromNCBI_path_to_bin_directory_CreateBD_MOD = TextFromNCBI_path_to_bin_directory + doubleBackSlash
            Script3.Path_To_NCBI_BLAST_Bin_Directory = TextFromNCBI_path_to_bin_directory_CreateBD_MOD

            Script3.Path_To_Database_Fasta_File = TextFrom_Path_to_Genome_fasta_file
            Script3.Data_Base_Type = TextFromDropDownMenu_DataBaseType

            print(Script3.Path_To_NCBI_BLAST_Bin_Directory)
            print(Script3.Path_To_Database_Fasta_File)
            print(Script3.Data_Base_Type)

            # Output Messages

            self.message = Label(text="Running", fg='red').grid(row=19, column=2, columnspan=2, sticky=W + E, padx=5,
                                                                pady=5)
            Script3.main(self)
            self.message1 = Label(text="Completed!", fg='red').grid(row=19, column=2, columnspan=2, sticky=W + E,
                                                                    padx=5,
                                                                    pady=5)
            messagebox.showinfo('Information', "Completed!")

        else:
            messagebox.showerror("Error", "All required parameters must be filled")
        self.wind.lift()

    def Select_button_for_Create_Database_for_NCBI_BLAST(self):
        Input_file_path_NCBI_Bin_Directory = filedialog.askdirectory(initialdir='.')
        print(Input_file_path_NCBI_Bin_Directory)
        self.NCBIPath_to_BIN.delete(0, END)
        self.NCBIPath_to_BIN.insert(0, Input_file_path_NCBI_Bin_Directory)
        self.wind.lift()
        return

    def select_path_for_FASTA_file(self):
        folder_selected_for_FASTA_File = filedialog.askopenfilename(initialdir='/', title="Select file",
                                                                    filetypes=(
                                                                        ("FASTA files", "*.fasta"),
                                                                        ("all files", "*.*")))
        print(folder_selected_for_FASTA_File)
        self.Path_to_FASTA_file.delete(0, END)
        self.Path_to_FASTA_file.insert(0, folder_selected_for_FASTA_File)
        self.wind.lift()
        return

    def Drop_down_definition_selected(self):
        return

    def ClearAll_Create_GENOME_DB(self):
        self.NCBIPath_to_BIN.delete(0, END)
        self.Path_to_FASTA_file.delete(0, END)

    def close_window(self):
        self.wind.destroy()


# BLAST and Filtering
class Win5(Win1):
    def __init__(self, window):
        # Initializations
        self.window = window
        self.wind = window
        self.wind.title("BLAST and filtering")
        self.wind.wm_iconbitmap(
            'Currito.ico')

        # Creating a Frame Container
        frame = LabelFrame(self.wind, text="Complete the following parameters")
        frame.grid(row=0, column=0, columnspan=3, padx=20, pady=20)

        # Select your NCBI bin Path in BF
        # Row 1
        Label(frame, text="Select path to NCBI/bin directory : ").grid(row=1, column=0, sticky=W)
        Label(frame, text="Example: " + "C:" + doubleBackSlash +
                          doubleBackSlash + "NCBI" + doubleBackSlash + "blast-2.8.0+" + doubleBackSlash + "bin   ").grid(
            row=1, column=6, sticky=W)
        self.NCBIPath_to_BIN_in_BF = tk.Entry(frame)
        self.NCBIPath_to_BIN_in_BF.focus()
        self.NCBIPath_to_BIN_in_BF.grid(row=1, column=1, columnspan=2, ipadx=200, padx=5, pady=5, )
        # Button Select
        tk.Button(frame, text="Select", command=self.Select_button_for_Select_NCBI_path_in_BF).grid(row=1,
                                                                                                    column=5,
                                                                                                    padx=2,
                                                                                                    pady=2)
        # Select your NCBI Path to database
        # Row 2 to row 4
        Label(frame, text="Select path to your NCBI database: ").grid(row=2, column=0, sticky=W)
        Label(frame, text="Example: " + "C:" + doubleBackSlash + doubleBackSlash
                          + "Users/MyDocuments/Genome/Genome.fasta ").grid(row=2, column=6, sticky=W)
        Label(frame, text="Note: There should be other files created when the database was made.").grid(row=3, column=6,
                                                                                                        sticky=W)
        Label(frame, text="Example: Genome.fasta.nhr / Genome.fasta.nin / Genome.fasta.nsq ").grid(row=4, column=6,
                                                                                                   sticky=W)

        self.SelectDataBase_BF = tk.Entry(frame)
        self.SelectDataBase_BF.focus()
        self.SelectDataBase_BF.grid(row=2, column=1, columnspan=2, ipadx=200, padx=5, pady=5, )
        # Button Select
        tk.Button(frame, text="Select", command=self.Select_button_for_Select_NCBI_Database_in_BF).grid(row=2, column=5,
                                                                                                        padx=2, pady=2)

        # Select path to the query
        # Row 5
        Label(frame, text="Select path to your query: ").grid(row=5, column=0, sticky=W)
        Label(frame, text="Example: " + "C:" + doubleBackSlash + doubleBackSlash
                          + "Users/MyDocuments/Genome/MyFile.fasta ").grid(row=5, column=6, sticky=W)
        self.SelectPathToQuery_in_BF = tk.Entry(frame)
        self.SelectPathToQuery_in_BF.grid(row=5, column=1, columnspan=2, ipadx=200, padx=5, pady=5, )
        # Button Select
        tk.Button(frame, text="Select", command=self.Select_button_for_Select_Query_in_BF).grid(row=5, column=5, padx=2,
                                                                                                pady=2)

        # Select output path in BF
        # Row 6
        Label(frame, text="Output path directory: ").grid(row=6, column=0, sticky=W)
        Label(frame, text="Example: " + "C:" + doubleBackSlash +
                          doubleBackSlash + "Users/MyDocuments/ ").grid(row=6, column=6, sticky=W)
        self.OutputPathFor_BF = tk.Entry(frame)
        self.OutputPathFor_BF.grid(row=6, column=1, columnspan=2, ipadx=200, padx=5, pady=5, )
        # Button Select
        tk.Button(frame, text="Select", command=self.Button_select_directory_output_path_BF).grid(row=6,
                                                                                                  column=5,
                                                                                                  padx=2, pady=2)

        # Select output name in BF
        # Row 7
        Label(frame, text="Output file name: ").grid(row=7, column=0, sticky=W)
        Label(frame, text=" Example: My_BLAST_results").grid(row=7, column=6, sticky=W)
        Label(frame, text="").grid(row=7, column=6, sticky=W)
        self.OutputFileName_BF = tk.Entry(frame)
        self.OutputFileName_BF.grid(row=7, column=1, columnspan=2, ipadx=200, padx=5, pady=5, )

        # Row 8
        # SELECT THE PATH TO THE FILE InputFile_with_unique_ID.txt
        # Row 8
        Label(frame, text="Path to file with Unique ID (Optional): ").grid(row=8, column=0, sticky=W)
        Label(frame, text="Example: My_BLAST_results_InputFile_with_unique_ID.txt ").grid(row=8, column=6, sticky=W)
        Label(frame, text=" ").grid(row=13, column=6, sticky=W)
        self.InputFile_UniqueID_BF = tk.Entry(frame)
        self.InputFile_UniqueID_BF.grid(row=8, column=1, columnspan=2, ipadx=200, padx=5, pady=5, sticky=W)
        # Button Select
        tk.Button(frame, text="Select", command=self.Button_select_file_Unique_ID).grid(row=8,
                                                                                        column=5,
                                                                                        padx=2, pady=2,
                                                                                        sticky=W)

        # Word Size
        # Row 9
        Label(frame, text="Word size: ").grid(row=9, column=0, sticky=W)
        Label(frame, text="Recommended value:  11 ").grid(row=9, column=2, sticky=W)
        self.WordSize_BF = tk.Entry(frame)
        self.WordSize_BF.grid(row=9, column=1, columnspan=1, ipadx=100, padx=5, pady=5, sticky=W)

        # SELECT YOUR PERCENTAGE IDENTITY = DEFAULT IS 70 %, you can change it. - f 70
        # Row 10
        Label(frame, text="Percentage Identity: ").grid(row=10, column=0, sticky=W)
        Label(frame, text="Recommended value:  70 ").grid(row=10, column=2, sticky=W)
        # Label(frame, text=" ").grid(row=9, column=6, sticky=W)
        self.PercentageIdentity_BF = tk.Entry(frame)
        self.PercentageIdentity_BF.grid(row=10, column=1, columnspan=1, ipadx=100, padx=5, pady=5, sticky=W)

        # SELECT YOUR THREADS. DEFAULT IS 4 IF YOUR COMPUTER HAS ONLY 4 CPU's, YOU CAN CHANGE it              -g 4
        # Row 11
        Label(frame, text="Number of threads: ").grid(row=11, column=0, sticky=W)
        Label(frame, text="Example: 4 ").grid(row=11, column=2, sticky=W)
        self.NumberOfThreads_BF = tk.Entry(frame)
        self.NumberOfThreads_BF.grid(row=11, column=1, columnspan=1, ipadx=100, padx=5, pady=5, sticky=W)

        # SELECT OUTPUT FORMAT. DEFAULT IS 6 AS A TABLE. YOU CAN SELECT OTHER FORMATS.                        -i 6
        # Row 12
        Label(frame, text="Output format: ").grid(row=12, column=0, sticky=W)
        Label(frame, text="Recommended format: 6 ").grid(row=12, column=2, sticky=W)
        self.OutputBLAST_Format_BF = tk.Entry(frame)
        self.OutputBLAST_Format_BF.grid(row=12, column=1, columnspan=1, ipadx=100, padx=5, pady=5, sticky=W)

        # SELECT YOUR PERCENTAGE OVERLAP. DEFAULT IS 80%.      -j 0.8
        # Row 13
        Label(frame, text="Percentage Overlap: ").grid(row=13, column=0, sticky=W)
        Label(frame, text="Recommended value: 0.8 ").grid(row=13, column=2, sticky=W)
        self.PercentageOverlap_Format_BF = tk.Entry(frame)
        self.PercentageOverlap_Format_BF.grid(row=13, column=1, columnspan=1, ipadx=100, padx=5, pady=5, sticky=W)

        # SELECT YOUR BITSCORE VALUE, DEFAULT IS 50. IF YOU ARE UNSURE DON'T USE THIS PARAMETER              -k 50
        # Row 14
        Label(frame, text="Bitscore: ").grid(row=14, column=0, sticky=W)
        Label(frame, text="Recommended value: 50 ").grid(row=14, column=2, sticky=W)
        self.Bitscore_BF = tk.Entry(frame)
        self.Bitscore_BF.grid(row=14, column=1, columnspan=1, ipadx=100, padx=5,
                              pady=5, sticky=W)

        # tick menu
        # Database type
        # Row 15
        Label(frame, text="Discontiguous Mega BLAST: ").grid(row=15, column=0, sticky=W)
        self.checked = tk.IntVar()
        self.CheckBox_BF = tk.Checkbutton(frame, text="dc-megablast", variable=self.checked, onvalue=1, offvalue=0)
        self.CheckBox_BF.grid(row=15, column=1, columnspan=1, padx=5, pady=5, sticky=W)

        # Button clear
        tk.Button(frame, text="Clear all", command=lambda: [self.Button_clear_all_Blast_and_Filtering()]).grid(row=17,
                                                                                                               column=6,
                                                                                                               columnspan=1,
                                                                                                               padx=5,
                                                                                                               pady=5,
                                                                                                               sticky=W + E)
        # Button Run
        tk.Button(frame, text="Run", command=self.Button_run_BF).grid(row=18, column=6, columnspan=1, padx=5, pady=5,
                                                                      sticky=W + E)

        # BUTTON Close
        tk.Button(frame, text="Close", command=self.close_window).grid(row=19, column=6, columnspan=1, padx=5, pady=5,
                                                                       sticky=W + E)

    def Button_run_BF(self):
        TextFrom_Path_To_NCBI_Bin_BF = self.NCBIPath_to_BIN_in_BF.get()
        TextFrom_DC_MEGABLAST_BF = self.checked.get()
        TextFrom_Path_To_Database_BF = self.SelectDataBase_BF.get()
        Text_from_Path_to_Query_BF = self.SelectPathToQuery_in_BF.get()
        TextFrom_Path_Output_BF = self.OutputPathFor_BF.get()
        TextFrom_OutputFile_name_BF = self.OutputFileName_BF.get()
        Text_From_InputFile_UniqueID_BF = self.InputFile_UniqueID_BF.get()
        TextFromWordSize_BF = self.WordSize_BF.get()
        TextFromPercentageID_BF = self.PercentageIdentity_BF.get()
        TextFromNumberOfThreads_BF = self.NumberOfThreads_BF.get()
        TextFromOutputBLAST_BF = self.OutputBLAST_Format_BF.get()
        TextFromPercentageOverlap_BF = self.PercentageOverlap_Format_BF.get()
        TextFromBitscore_BF = self.Bitscore_BF.get()

        if (len(TextFrom_Path_To_NCBI_Bin_BF) != 0

                and len(TextFrom_Path_To_Database_BF) != 0
                and len(Text_from_Path_to_Query_BF) != 0
                and len(TextFrom_Path_Output_BF) != 0
                and len(TextFrom_OutputFile_name_BF) != 0
                and len(TextFromWordSize_BF) != 0
                and len(TextFromPercentageID_BF) != 0
                and len(TextFromNumberOfThreads_BF) != 0
                and len(TextFromOutputBLAST_BF) != 0
                and len(TextFromPercentageOverlap_BF) != 0
                and len(TextFromBitscore_BF) != 0):

            # Here
            TextFrom_Path_To_NCBI_Bin_BF_MOD = TextFrom_Path_To_NCBI_Bin_BF + doubleBackSlash
            Script4.x_Path_to_NCBI_Directory_BF = TextFrom_Path_To_NCBI_Bin_BF_MOD

            Script4.y_DC_MegaBlast_BF = TextFrom_DC_MEGABLAST_BF

            TextFrom_Path_To_Database_BF_MOD = TextFrom_Path_To_Database_BF #+ doubleBackSlash
            Script4.a_Data_Base_fasta = TextFrom_Path_To_Database_BF_MOD

            Text_from_Path_to_Query_BF_MOD = Text_from_Path_to_Query_BF #+ doubleBackSlash
            Script4.b_Query_fasta_file = Text_from_Path_to_Query_BF_MOD

            TextFrom_Path_Output_BF_MOD = TextFrom_Path_Output_BF + doubleBackSlash
            Script4.c_Output_Path_ = TextFrom_Path_Output_BF_MOD

            Script4.d_Output_file_name = TextFrom_OutputFile_name_BF
            Script4.e_word_size = TextFromWordSize_BF
            Script4.f_Percentage_identity = TextFromPercentageID_BF
            Script4.g_number_of_threads = TextFromNumberOfThreads_BF
            Script4.i_OutputFormat = TextFromOutputBLAST_BF
            Script4.j_Percentage_overlap = TextFromPercentageOverlap_BF
            Script4.k_bitscore = TextFromBitscore_BF
            Script4.l_InputFile_with_unique_ID = Text_From_InputFile_UniqueID_BF

            print(Script4.x_Path_to_NCBI_Directory_BF)
            print("")
            print(Script4.a_Data_Base_fasta)
            print(Script4.b_Query_fasta_file)
            print(Script4.c_Output_Path_)
            print(Script4.d_Output_file_name)
            print(Script4.e_word_size)
            print(Script4.f_Percentage_identity)
            print(Script4.g_number_of_threads)
            print(Script4.i_OutputFormat)
            print(Script4.j_Percentage_overlap)
            print(Script4.k_bitscore)
            print(Script4.l_InputFile_with_unique_ID)
            print(Script4.y_DC_MegaBlast_BF)

            # Output Messages

            self.message = Label(text="Running", fg='red').grid(row=19, column=2, columnspan=2, sticky=W + E, padx=5,
                                                                pady=5)
            Script4.main(self)
            self.message1 = Label(text="Completed!", fg='red').grid(row=19, column=2, columnspan=2, sticky=W + E,
                                                                    padx=5,
                                                                    pady=5)
            messagebox.showinfo('Information', "Completed!")
        else:
            messagebox.showerror("Error", "All required parameters must be filled")
        self.wind.lift()

    def Select_button_for_Select_NCBI_path_in_BF(self):
        Input_file_path_NCBI_Bin_Directory_BF = filedialog.askdirectory(initialdir='.')
        print(Input_file_path_NCBI_Bin_Directory_BF)
        self.NCBIPath_to_BIN_in_BF.delete(0, END)
        self.NCBIPath_to_BIN_in_BF.insert(0, Input_file_path_NCBI_Bin_Directory_BF)
        self.wind.lift()
        return

    def Select_button_for_Select_NCBI_Database_in_BF(self):
        Input_file_Database_BF = filedialog.askopenfilename(initialdir='/', title="Select file",
                                                            filetypes=(
                                                                ("Fasta files", "*.fasta"), ("all files", "*.*")))
        print(Input_file_Database_BF)
        self.SelectDataBase_BF.delete(0, END)
        self.SelectDataBase_BF.insert(0, Input_file_Database_BF)
        self.wind.lift()
        return

    def Select_button_for_Select_Query_in_BF(self):
        Input_file_query_BF = filedialog.askopenfilename(initialdir='/', title="Select file",
                                                         filetypes=(("Fasta files", "*.fasta"), ("all files", "*.*")))
        print(Input_file_query_BF)
        self.SelectPathToQuery_in_BF.delete(0, END)
        self.SelectPathToQuery_in_BF.insert(0, Input_file_query_BF)
        self.wind.lift()
        return

    def Button_select_directory_output_path_BF(self):
        Output_path_in_BF = filedialog.askdirectory(initialdir='.')
        print(Output_path_in_BF)
        self.OutputPathFor_BF.delete(0, END)
        self.OutputPathFor_BF.insert(0, Output_path_in_BF)
        self.wind.lift()
        return

    def Button_select_file_Unique_ID(self):
        Path_to_file_with_Unique_ID = filedialog.askopenfilename(initialdir='/', title="Select file", filetypes=(
            ("Text files", "*.txt"), ("all files", "*.*")))
        print(Path_to_file_with_Unique_ID)
        self.InputFile_UniqueID_BF.delete(0, END)
        self.InputFile_UniqueID_BF.insert(0, Path_to_file_with_Unique_ID)
        self.wind.lift()
        return

    def Button_clear_all_Blast_and_Filtering(self):
        self.NCBIPath_to_BIN_in_BF.delete(0, END)
        self.SelectDataBase_BF.delete(0, END)
        self.SelectPathToQuery_in_BF.delete(0, END)
        self.OutputPathFor_BF.delete(0, END)
        self.OutputFileName_BF.delete(0, END)
        self.InputFile_UniqueID_BF.delete(0, END)
        self.PercentageIdentity_BF.delete(0, END)
        self.NumberOfThreads_BF.delete(0, END)
        self.OutputBLAST_Format_BF.delete(0, END)
        self.PercentageOverlap_Format_BF.delete(0, END)
        self.Bitscore_BF.delete(0, END)
        self.WordSize_BF.delete(0, END)

    def close_window(self):
        self.wind.destroy()









#Filtering
class Win55(Win1):
    def __init__(self, window):
        # Initializations
        self.window = window
        self.wind = window
        self.wind.title("Filtering")
        self.wind.wm_iconbitmap(
            'Currito.ico')

        # Creating a Frame Container
        frame = LabelFrame(self.wind, text="Complete the following parameters")
        frame.grid(row=0, column=0, columnspan=3, padx=20, pady=20)

 # Select path to the query in F
        # Row 2
        Label(frame, text="Select path to your BLAST file(*): ").grid(row=2, column=0, sticky=W)
        Label(frame, text="Example: " + "C:" + doubleBackSlash + doubleBackSlash
                          + "Users/MyDocuments/MyBLAST_results.txt ").grid(row=2, column=6, sticky=W)
        self.SelectPathToQuery_in_F = tk.Entry(frame)
        self.SelectPathToQuery_in_F.grid(row=2, column=1, columnspan=2, ipadx=200, padx=5, pady=5, )
        # Button Select
        tk.Button(frame, text="Select", command=self.Select_button_for_Select_Query_in_F).grid(row=2, column=5, padx=2,
                                                                                                pady=2)

        # Select output path in F
        # Row 3
        Label(frame, text="Output path directory: ").grid(row=3, column=0, sticky=W)
        Label(frame, text="Example: " + "C:" + doubleBackSlash +
                          doubleBackSlash + "Users/MyDocuments/ ").grid(row=3, column=6, sticky=W)
        self.OutputPathFor_F = tk.Entry(frame)
        self.OutputPathFor_F.grid(row=3, column=1, columnspan=2, ipadx=200, padx=5, pady=5, )
        # Button Select
        tk.Button(frame, text="Select", command=self.Button_select_directory_output_path_F).grid(row=3,
                                                                                                  column=5,
                                                                                                  padx=2, pady=2)

        # Select output name in F
        # Row 4
        Label(frame, text="Output file name: ").grid(row=4, column=0, sticky=W)
        Label(frame, text=" Example: My_filtered_BLAST_results").grid(row=4, column=6, sticky=W)
        Label(frame, text="").grid(row=4, column=6, sticky=W)
        self.OutputFileName_F = tk.Entry(frame)
        self.OutputFileName_F.grid(row=4, column=1, columnspan=2, ipadx=200, padx=5, pady=5, )




        # Row 6
        Label(frame, text="Percentage Overlap: ").grid(row=6, column=0, sticky=W)
        Label(frame, text="Recommended value: 0.8 ").grid(row=6, column=2, sticky=W)
        self.PercentageOverlap_Format_F = tk.Entry(frame)
        self.PercentageOverlap_Format_F.grid(row=6, column=1, columnspan=1, ipadx=100, padx=5, pady=5, sticky=W)

        # SELECT YOUR BITSCORE VALUE, DEFAULT IS 50. IF YOU ARE UNSURE DON'T USE THIS PARAMETER              -k 50
        # Row 7
        Label(frame, text="Bitscore: ").grid(row=7, column=0, sticky=W)
        Label(frame, text="Recommended value: 50 ").grid(row=7, column=2, sticky=W)
        self.Bitscore_F = tk.Entry(frame)
        self.Bitscore_F.grid(row=7, column=1, columnspan=1, ipadx=100, padx=5,
                              pady=5, sticky=W)

        #Rown 8 Notes
        Label(frame, text="").grid(row=8, column=0, sticky=W)
        Label(frame, text="* The BLASTn output format: ").grid(row=9, column=0, sticky=W)
        Label(frame, text=" TABULAR OUTPUT FORMAT: 6").grid(row=10, column=0, sticky=W)
        Label(frame, text=" COLUMN HEADERS:").grid(row=11, column=0, sticky=W)
        Label(frame, text="qseqid sacc stitle qseq sseq ").grid(row=12, column=0, sticky=W)
        Label(frame, text="nident mismatch pident length ").grid(row=13, column=0, sticky=W)
        Label(frame, text=" evalue bitscore qstart qend sstart send").grid(row=14, column=0, sticky=W)
        Label(frame, text=" gapopen gaps qlen slen").grid(row=15, column=0, sticky=W)
        # Button clear
        tk.Button(frame, text="Clear all", command=lambda: [self.Button_clear_all_Filtering()]).grid(row=17,
                                                                                                           column=6,
                                                                                                           columnspan=1,
                                                                                                           padx=5,
                                                                                                           pady=5,
                                                                                                           sticky=W + E)
        # Button Run
        tk.Button(frame, text="Run", command=self.Button_run_F).grid(row=18, column=6, columnspan=1, padx=5, pady=5,
                                                                  sticky=W + E)

        # BUTTON Close
        tk.Button(frame, text="Close", command=self.close_window).grid(row=19, column=6, columnspan=1, padx=5, pady=5,
                                                                   sticky=W + E)


    def Button_run_F(self):
        Text_from_Path_to_Query_F = self.SelectPathToQuery_in_F.get()
        TextFrom_Path_Output_F = self.OutputPathFor_F.get()
        TextFrom_OutputFile_name_F = self.OutputFileName_F.get()
        #TextFromPercentageID_F = self.PercentageIdentity_F.get()
        TextFromPercentageOverlap_F = self.PercentageOverlap_Format_F.get()
        TextFromBitscore_F = self.Bitscore_F.get()

        if (len(Text_from_Path_to_Query_F) != 0

                and len(TextFrom_Path_Output_F) != 0
                and len(TextFrom_OutputFile_name_F) != 0
                and len(TextFromPercentageOverlap_F) != 0
                and len(TextFromBitscore_F) != 0):

            # Here

            Text_from_Path_to_Query_F_MOD = Text_from_Path_to_Query_F #+ doubleBackSlash
            Script55.a_BLAST_input_path_and_file_ = Text_from_Path_to_Query_F_MOD

            TextFrom_Path_Output_F_MOD = TextFrom_Path_Output_F + doubleBackSlash
            Script55.b_Output_Path_ = TextFrom_Path_Output_F_MOD

            Script55.c_Output_file_name = TextFrom_OutputFile_name_F
            #Script55.f_Percentage_identity = TextFromPercentageID_F
            Script55.d_Percentage_overlap = TextFromPercentageOverlap_F
            Script55.e_bitscore = TextFromBitscore_F



            print(Script55.a_BLAST_input_path_and_file_)
            print(Script55.b_Output_Path_)
            print(Script55.c_Output_file_name)
            #print(Script55.f_Percentage_identity)
            print(Script55.d_Percentage_overlap)
            print(Script55.e_bitscore)


            # Output Messages

            self.message = Label(text="Running", fg='red').grid(row=19, column=2, columnspan=2, sticky=W + E, padx=5,
                                                                pady=5)
            Script55.main(self)
            self.message1 = Label(text="Completed!", fg='red').grid(row=19, column=2, columnspan=2, sticky=W + E,
                                                                    padx=5,
                                                                    pady=5)
            messagebox.showinfo('Information', "Completed!")
        else:
            messagebox.showerror("Error", "All required parameters must be filled")
        self.wind.lift()


    def Select_button_for_Select_Query_in_F(self):
        Input_file_query_F = filedialog.askopenfilename(initialdir='/', title="Select file",
                                                         filetypes=(("Text file", "*.txt"), ("all files", "*.*")))
        print(Input_file_query_F)
        self.SelectPathToQuery_in_F.delete(0, END)
        self.SelectPathToQuery_in_F.insert(0, Input_file_query_F)
        self.wind.lift()
        return

    def Button_select_directory_output_path_F(self):
        Output_path_in_F = filedialog.askdirectory(initialdir='.')
        print(Output_path_in_F)
        self.OutputPathFor_F.delete(0, END)
        self.OutputPathFor_F.insert(0, Output_path_in_F)
        self.wind.lift()
        return

    def Button_clear_all_Filtering(self):

        self.SelectPathToQuery_in_F.delete(0, END)
        self.OutputPathFor_F.delete(0, END)
        self.OutputFileName_F.delete(0, END)
        self.PercentageOverlap_Format_F.delete(0, END)
        self.Bitscore_F.delete(0, END)


    def close_window(self):
        self.wind.destroy()











# Help or instructions
class Win6(Win1):
    def __init__(self, window):
        # Initializations
        self.window = window
        self.wind = window
        self.wind.title("Instructions")
        self.wind.wm_iconbitmap(
            'Currito.ico')

        # Creating a Frame Container
        frame = LabelFrame(self.wind, text="User guide")
        frame.grid(row=0, column=0, columnspan=3, padx=20, pady=20)

        # Select your NCBI Path
        Label(frame, text="").grid(row=1, column=0)

        Label(frame,
              text="Installation requirements:").grid(
            row=2, column=0, sticky=W)
        #Label(frame, text="    *Python 3 or above (https://www.python.org/downloads/)").grid(row=3, column=0, sticky=W)
        #3Label(frame, text="    *Python module BioPython (https://biopython.org/wiki/Download) ").grid(row=4, column=0, sticky=W)
        Label(frame,
              text="    *BLAST+ greater than v2.6 or the latest version (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)").grid(
            row=5, column=0, sticky=W)

        Label(frame,
              text="").grid(
            row=6, column=0, sticky=W)

        Label(frame,
              text="Recommendations:").grid(
            row=8, column=0, sticky=W)

        #Label(frame, text="    *Add Python to the environment path to be able to run this script.").grid( row=9, column=0, sticky=W)
        #Label(frame, text="        -To add python to the path in Windows, you can do it by modifying it in: ").grid(row=10, column=0, sticky=W)
        #Label(frame, text="        Control Panel > System and Security > System > Advanced System Settings > Environment Variables > System Variables > Path").grid( row=11, column=0, sticky=W)
        Label(frame,
              text="    *Do not install BLAST in the 'Program Files' directory. The space between words will make this script to crash.").grid(
            row=12, column=0, sticky=W)
        Label(frame,
              text="    *Do not save any document using names with spaces between words, use underscores. Example: My_file.").grid(
            row=13, column=0, sticky=W)

        Label(frame, text="    *The NCBI BLAST also takes .fna files, in addition to .fasta files.").grid(
            row=14, column=0, sticky=W)

        Label(frame,
              text=" ").grid(
            row=15, column=0, sticky=W)
        Label(frame,
              text="").grid(
            row=16, column=0, sticky=W)

        Label(frame, text="Questions or comments: berenicetalamantes@yahoo.fr").grid(
            row=20, column=0, sticky=W)

        Label(frame, text="Developed by : Berenice Talamantes-Becerra, Jason Carling, Arthur Georges").grid(
            row=21, column=0, sticky=W)

        # BUTTON Close
        tk.Button(frame, text="Close", command=self.close_window).grid(row=23, column=2, columnspan=1, padx=5, pady=5,
                                                                       sticky=W + E)

    def close_window(self):
        self.wind.destroy()


Script_1_Instance = Script1()
if __name__ == "__main__":
    window = Tk()
    application = Win1(window)
    window.mainloop()
