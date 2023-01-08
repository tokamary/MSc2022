# 
#                                    ORF DNA Sequence Finder (ODSF)

# This is a python 3 programm. It takes a multifasta file or a sequence-string and returns, if there are, coding regions 
# (ORF) and its length.
# First, the programm ask you to enter if you want to read from a multifasta file (by writing 'file') or to type one 
# sequence (by writing 'seq').
# After that, based on your choice, you need to type the name of the file or your sequence.
# The program  can reject DNA sequences with degenerative bases and can handle 
# input with spaces,'-' as well as upper/lower case letters. 
# To excecute it from command line you have to install python 3 and give the commant: python ODSF.py
##------------------------------------------------------------------------------------------------------------
# Submitted by Tokamani Maria
# MSc 'Applied Bioinformatics and Data Analysis' DUTh, Alexandroupolis, GR.
# 2022-2023
# Module EBP101 - Unix and Programming Basics, N.M.Glykos
##------------------------------------------------------------------------------------------------------------

import re       #"re" module used for string searching and manipulation
import string   # handling strings
import sys      #


# function to read multiple fasta files
def read_fasta(fp):
        name = None
        seq = []
        for line in fp:
            line = line.rstrip()
            if line.startswith(">"):
                if name: yield (name, ''.join(seq))
                name = line
                seq = []
            else:
                seq.append(line)
        if name: yield (name, ''.join(seq))

# Below the program uses the test sequence to run directly. 
#seq_1=('ttaattaCGTGCTGCTAtG  GCTaTAaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaTGCTGTGCTGACatgGCTTACTCTGTAaA  CTATtTACAtgc')
#seq_1=('MALwMRLLPL LALLALWGPD PAAaFVNQHL CGSHLVEALY LVCGERGFFY TPKTRREAED LQVGqVELGG GPgAGSLQPL ALEGsLQKRG IVQCCTSIC SLYQLENYCN') #checked and it works!

#solution
#number strand frame start end len [>30bp]
#ORF1   +   1   X
#ORF2	+	2	17	115	99
#ORF3   +   3           12
#ORF4	-	1	101	21	81
#ORF5   -   2   X
#ORF6   -   3   1   4   6

# User can input a sequence from a file. To check this, comment out the lines below


UserChoice = input('Do you want the programm read sequence from a file or to type your DNA sequence on the terminal. [For file type "file", for a sequence type "seq"]: ') 

if UserChoice == 'file':
    filename = input('The programm operates only with fasta format sequence files with one or multiple sequences. Copy the fasta file into the same folder with the programm. Type here the name of the file: ')
  
    sequence = []
    	  
    with open(filename) as fp:
        for name, seq in read_fasta(fp):
            sequence.append((name, seq))
        #print(sequence)

            for i in range (0,len(sequence),1):
                allseq= []
                allname=[]
                allseq.append(sequence[i][1]) #retrieve only sequences
                allname.append(sequence[i][0]) #retrieve only the name description

            #print(allname)


            for x in range(0,len(allseq),1):
                seq_1 = allseq[x] #make a string one sequences
                name = allname[x]
                listOFDbases=  re.findall("[^AaTtGgCc\s-]", seq_1) #detect aminoacids or degenerate bases
                numberofDbases= len(listOFDbases)

                if numberofDbases > 0:
                    print ("@ - Sequence under evaluation: ", name, ". Degenerate bases or aminoacid sequence was detected!")
                    
                   

                else:
                    print ("@ -> A DNA sequence was detected. Searching for possible ORF...Please wait.")

                    # reads, 'cleans' and prepare the sequence.
                    capital_seq = seq_1.upper() # converts all to upper
                    profinalseq = capital_seq.replace("-","") # eliminates '-'
                    finalseq = profinalseq.replace(" ", "") # eliminates gaps
                    
                    #creates the cDNA strand with 5'-3'.
                    
                    cDNA_seq = finalseq.translate(str.maketrans('ATGC', 'TACG'))  #replace nucleotides based on complimentarity rule
                    final_cDNA = cDNA_seq [::-1] #reverse the orientation of the sequence for 3-5 to 5-3

                    print ("The sequence is: ", name )
                    print ("The 5' - 3' DNA sequence of your input is: \n", finalseq )
                    print ("The 5' - 3' cDNA sequence of your input is: \n", final_cDNA )
                    print ("The total length of the sequence is:", len(finalseq))

                    # create the 3 frames of the sequence + strand
                    frames = []  #list that include the 3 + frames
                    frames_rev = [] #list that include the 3 - frames

                    # take into account the 3 different frames in each strand
                    
                    for y in range(0,3):

                        #tmp list of frames for an open frame  
                        frame_tmp=[]  #frames of the + strand
                        frame_tmp_rev=[]  #frames of the - strand

                        #split the sequence string on triplets, and append to list

                        for x in range(y,len(finalseq)+1,3):

                            if x+2 <= len(finalseq)-1:

                                tmp_tri= finalseq[x]+ finalseq[x+1] + finalseq[x+2]        

                                frame_tmp.append(tmp_tri)
                        frames.append(frame_tmp)

                        for z in range(y,len(final_cDNA)+1,3):

                            if z+2 <= len(final_cDNA)-1:

                                tmp_tri_rev= final_cDNA[z]+ final_cDNA[z+1] + final_cDNA[z+2]        

                                frame_tmp_rev.append(tmp_tri_rev)
                        frames_rev.append(frame_tmp_rev)
                        
                    
                    #print(frames)
                    #print(frames_rev)
                    
                
                    #Find ATG and stop condon and extract the ORF
                    listOfOrf=[]
                    listOfOrf_rev=[]
                    no = 1

                    for i in range(0,len(frames),1): #looping all the frames of + strand
                        start=0
                        while start < len(frames[i]): #looping each frame for start and stop codons 
                            if frames[i][start]=="ATG":
                                for stop in range(start+1,len(frames[i]),1):
                                            if frames[i][stop]=="TAA" or  frames[i][stop]=="TAG" or  frames[i][stop]=="TGA" :
                                                listOfOrf.append(frames[i][start:(stop+1)]) # retrieve the orf 
                                                orf_list = (frames[i][start:(stop+1)]) # one line list
                                                orf_str = ' '.join(orf_list)  #convert list to string
                                                orf_len = len(orf_list)*3 #calculate length of the coding region
                                                print(">ORF was detected on the + strand. Number", (no))
                                                print("Its lenght is:", orf_len )
                                                print("and the sequence is: \n", orf_str )
                                                no = no +1
                                                start=stop+1 # avoiding multiple start codons
                                                break
                            start+=1

                    
                    for i in range(0,len(frames_rev),1): #looping all the frames of - strand
                        start=0
                        while start <len(frames_rev[i]): #looping each frame for start and stop codons 
                            if frames_rev[i][start]=="ATG":
                                for stop in range(start+1,len(frames_rev[i]),1):
                                            if frames_rev[i][stop]=="TAA" or  frames_rev[i][stop]=="TAG" or  frames_rev[i][stop]=="TGA" :
                                                listOfOrf_rev.append(frames_rev[i][start:(stop+1)]) # retrieve the orf 
                                                orf_list = (frames_rev[i][start:(stop+1)])
                                                orf_str = ' '.join(orf_list)
                                                orf_len = len(orf_list)*3
                                                print(">ORF was detected on the - strand. Number", (no))
                                                print("Its lenght is:", orf_len )
                                                print("and the sequence is: \n", orf_str )
                                                no = no +1
                                                start=stop+1 # avoiding multiple start codons
                                                break
                            start+=1
                
                    if len(listOfOrf) == 0:
                        print (">Zero ORF was found in the + strand of your DNA sequence")
                    
                    if len(listOfOrf_rev) == 0:
                        print (">Zero ORF was found in the - strand of your DNA sequence")
            

if UserChoice =='seq':


    all_seq = []

    seq_1 = ''

    #start a loop that will run until the user exit the programm

    while seq_1 != 'exit' :

        seq_1 = input('Type a DNA sequence with orientation 5΄ - 3΄ and press enter once. [If you want to exit the programm type "exit" and press enter.]: ') 

        if seq_1 != 'exit':

            # reject the sequence if detects degenerative bases or aminoacid sequence.
            listOFDbases=  re.findall("[^AaTtGgCc\s-]", seq_1)
            numberofDbases= len(listOFDbases)

            if numberofDbases > 0:
                print("Degenerate bases or aminoacid sequence was detected!")  
                

            else:
                print ("@ -> A DNA sequence was detected. Searching for possible ORF...Please wait.")

                # reads, 'cleans' and prepare the sequence.
                capital_seq = seq_1.upper() # converts all to upper
                profinalseq = capital_seq.replace("-","") # eliminates '-'
                finalseq = profinalseq.replace(" ", "") # eliminates gaps

                #creates the cDNA strand with 5'-3'.
                
                cDNA_seq = finalseq.translate(str.maketrans('ATGC', 'TACG'))  #replace nucleotides based on complimentarity rule
                final_cDNA = cDNA_seq [::-1] #reverse the orientation of the sequence for 3-5 to 5-3

                print ("The 5' - 3' DNA sequence of your input is: \n", finalseq )
                print ("The 5' - 3' cDNA sequence of your input is: \n", final_cDNA )
                print ("The total length of the sequence is:", len(finalseq))

                # create the 3 frames of the sequence + strand
                frames = []  #list that include the 3 + frames
                frames_rev = [] #list that include the 3 - frames

                # take into account the 3 different frames in each strand
                
                for y in range(0,3):

                    #tmp list of frames for an open frame  
                    frame_tmp=[]  #frames of the + strand
                    frame_tmp_rev=[]  #frames of the - strand

                    #split the sequence string on triplets, and append to list

                    for x in range(y,len(finalseq)+1,3):

                        if x+2 <= len(finalseq)-1:

                            tmp_tri= finalseq[x]+ finalseq[x+1] + finalseq[x+2]        

                            frame_tmp.append(tmp_tri)
                    frames.append(frame_tmp)

                    for z in range(y,len(final_cDNA)+1,3):

                        if z+2 <= len(final_cDNA)-1:

                            tmp_tri_rev= final_cDNA[z]+ final_cDNA[z+1] + final_cDNA[z+2]        

                            frame_tmp_rev.append(tmp_tri_rev)
                    frames_rev.append(frame_tmp_rev)
                    
                
                #print(frames)
                #print(frames_rev)
                
            
                #Find ATG and stop condon and extract the ORF
                listOfOrf=[]
                listOfOrf_rev=[]
                no = 1

                for i in range(0,len(frames),1): #looping all the frames of + strand
                    start=0
                    while start < len(frames[i]): #looping each frame for start and stop codons 
                        if frames[i][start]=="ATG":
                            for stop in range(start+1,len(frames[i]),1):
                                        if frames[i][stop]=="TAA" or  frames[i][stop]=="TAG" or  frames[i][stop]=="TGA" :
                                            listOfOrf.append(frames[i][start:(stop+1)]) # retrieve the orf 
                                            orf_list = (frames[i][start:(stop+1)]) # one line list
                                            orf_str = ' '.join(orf_list)  #convert list to string
                                            orf_len = len(orf_list)*3 #calculate length of the coding region
                                            print(">ORF was detected on the + strand. Number", (no))
                                            print("Its lenght is:", orf_len )
                                            print("and the sequence is: \n", orf_str )
                                            no = no +1
                                            start=stop+1 # avoiding multiple start codons
                                            break
                        start+=1

                
                for i in range(0,len(frames_rev),1): #looping all the frames of - strand
                    start=0
                    while start <len(frames_rev[i]): #looping each frame for start and stop codons 
                        if frames_rev[i][start]=="ATG":
                            for stop in range(start+1,len(frames_rev[i]),1):
                                        if frames_rev[i][stop]=="TAA" or  frames_rev[i][stop]=="TAG" or  frames_rev[i][stop]=="TGA" :
                                            listOfOrf_rev.append(frames_rev[i][start:(stop+1)]) # retrieve the orf 
                                            orf_list = (frames_rev[i][start:(stop+1)])
                                            orf_str = ' '.join(orf_list)
                                            orf_len = len(orf_list)*3
                                            print(">ORF was detected on the - strand. Number", (no))
                                            print("Its lenght is:", orf_len )
                                            print("and the sequence is: \n", orf_str )
                                            no = no +1
                                            start=stop+1 # avoiding multiple start codons
                                            break
                        start+=1
            
                if len(listOfOrf) == 0:
                    print (">Zero ORF was found in the + strand of your DNA sequence")
                
                if len(listOfOrf_rev) == 0:
                    print (">Zero ORF was found in the - strand of your DNA sequence")
            
