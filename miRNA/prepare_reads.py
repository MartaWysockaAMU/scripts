import os
import sys
import gzip
import numpy as np
from Bio import SeqIO
from functools import partial
from mimetypes import guess_type
from prettytable import PrettyTable
from prettytable import PLAIN_COLUMNS
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqIO.QualityIO import FastqGeneralIterator


class Reads():

    def __init__(self, inp: str):
        """Constructor for class Reads

        Arguments:
            inp {str} -- reads input argument
            out {str} -- prefix for output file
            dir {str} -- output directory name
        """
        self.input = inp
        #self.out_pref = out
        self.samples = {}  # names of sample for paths
        self.reads_files = {}  # reads files splited by extension
		
    def check_format(self, file_path): 
        encoding = guess_type(file_path)[1]
        _open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open
        try:
            #if file is fastq
            with _open(file_path) as fq:
                record = next(SeqIO.parse(fq, "fastq"))
            return [True, "fastq"]
        except ValueError:
            try:
                #if file is fasta
                with _open(file_path) as handle:
                    record = next(SeqIO.parse(handle, "fasta"))
            except StopIteration:
                #if file is txt with paths
                line = _open(file_path).readline()
                if len(line.split()) > 1:
                    path = line.split()[0]
                    if os.path.exists(path):
                        return [True, "txt"]
                    else:
                        msg = "{arg} -> incorrect file type. Path is not exists\n".format(arg=self.input[0])
                        self.exception_handling(msg) 
                elif len(line.split()) == 1:
                    path = line.strip()
                    if os.path.exists(path):
                        return [True, "txt"]
                    else:
                        msg = "{arg} -> incorrect file type. Path is not exists\n".format(arg=self.input[0])
                        self.exception_handling(msg)
                else:
                    msg = "{arg} -> incorrect file type. Provide a fasta or fastq file for reads or txt file for paths to read files\n".format(arg=self.input[0])
                    self.exception_handling(msg)
            #if it is fasta check if it is collapsed
            try:
                int(record.name.strip().split("_x")[-1])  
                return [True, ("fasta", "_x")]
            except ValueError:
                try:
                    int(record.name.strip().split("-")[-1])
                    return [True, ("fasta", "-")] 
                except ValueError:
                    return [True, "fasta"]
           
    def get_and_split_paths(self):
        """Chcek if passed files exist and are not empty
           Get sample name for each path

        Raises:
            Exception: incorrect file type
        """
        self.reads_files = {"fasta": [], "fastq": [], ("fasta", "-"): [], ("fasta", "_x"): []}
        samples = {}  # names of sample for paths
        name = path = ""        
        if len(self.input) == 1:
            # if path to directory with files is provided
            if os.path.exists(self.input[0]):
                if os.path.isdir(self.input[0]):
                    path = os.path.split("{}/".format(self.input[0]))[0]
                    for read in os.listdir(path):
                        path_to_file = os.path.join(path, read)
                        if self.is_file(path_to_file) and not self.is_empty(path_to_file):
                            form = self.check_format(path_to_file)
                            if form[0]:
                                if '.gz' in read:
                                    name = ".".join(read.split(".")[0:-2])
                                else:
                                    name = ".".join(read.split(".")[0:-1])
                                samples[path_to_file] = name
                                self.reads_files[form[1]].append(path_to_file)
                elif self.is_file(self.input[0]) and\
                not self.is_empty(self.input[0]):
                    (path,name) = os.path.split("{}".format(self.input[0]))
                    path_to_file = os.path.join(path, name)
                    form = self.check_format(path_to_file)
                    
                    # if only one file with reads is provided
                    if form[0]:
                        if form[1] != "txt":
                            if '.gz' in name:
                                samples[path_to_file] = ".".join(name.split(".")[0:-2])
                            else:
                                samples[path_to_file] = ".".join(name.split(".")[0:-1])
                            self.reads_files[form[1]].append(path_to_file)
                        # if path to file with paths to reads file is provided
                        else:
                            encoding = guess_type(self.input[0])[1]  # uses file extension
                            _open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open
                            with _open(self.input[0]) as path_file:
                                for line in path_file:
                                    # if path and name for reads file is provided
                                    if len(line.split()) > 1:
                                        path = line.split()[0]
                                        if self.is_file(path) and\
                                        not self.is_empty(path):
                                            form = self.check_format(path)
                                            if form[0]:
                                                name = line.strip().split()[1]
                                                samples[path] = name
                                                self.reads_files[form[1]].append(path)
                                    # if only path is provided
                                    else:
                                        path_to_file = line.strip()
                                        name = os.path.split("{}".format(path_to_file))[1]
                                        if self.is_file(path_to_file) and\
                                        not self.is_empty(path_to_file):
                                            form = self.check_format(path_to_file)
                                            if form[0]:
                                                if '.gz' in name:
                                                    samples[path_to_file] = ".".join(name.split(".")[0:-2])
                                                else:
                                                    samples[path_to_file] = ".".join(name.split(".")[0:-1])
                                                self.reads_files[form[1]].append(path_to_file)
                else:
                    msg = "Provided path {p} doesn't exist\n".format(p=self.input[0])
                    self.exception_handling(msg) 
            else:
                msg = "Provided path {p} doesn't exist\n".format(p=self.input[0])
                self.exception_handling(msg)
        else:  # if reads files name are listed in arguments
            for path_to_file in self.input: 
                if self.is_file(path_to_file) and not self.is_empty(path_to_file):
                    form = self.check_format(path_to_file)
                    if form[0]:
                        name = name = os.path.split("{}".format(path_to_file))[1]
                        if '.gz' in name:
                            samples[path_to_file] = ".".join(name.split(".")[0:-2])
                        else:
                            samples[path_to_file] = ".".join(name.split(".")[0:-1])
                        self.reads_files[form[1]].append(path_to_file)
        self.samples = dict(samples)
        samples = {}
        
    def merge_reads(self, save: bool, out_file: str, save_stat: bool, gz: bool, plainTab: bool):
        """Collapse reads if not collapsed\
           and merges all collapsed reads files
        """
        sep = ""    
        if gz:    
            out_file = out_file + ".gz"
            merged = gzip.open(out_file, "wt")
        else:
            merged = open(out_file, "w")
        fq_count = len(self.reads_files["fastq"])
        coll_fa_count = len(self.reads_files[("fasta", "-")]) + len(self.reads_files[("fasta", "_x")])
        uncoll_fa_count = len(self.reads_files["fasta"])
        form_dict = {"Uncollapsed FASTA": uncoll_fa_count, "Collapsed FASTA": coll_fa_count, "FASTQ": fq_count}
        reads_count = {}
        for form in self.reads_files:
            for in_file in self.reads_files[form]:
                reads_count[in_file] = [] #[0] liczba unikalnych odczytów, [1] liczba wszystkich odczytów, [2] średnia długość liczona z wszsytkich (a może z unikalnych)
                uniq_reads = 0
                all_reads = 0
                lengths = []
                sample = self.samples[in_file]
                encoding = guess_type(in_file)[1]  # uses file extension
                _open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open
                if isinstance(form, tuple):
                    sep = form[1]
                    with _open(in_file) as fa:
                        for title, seq in SimpleFastaParser(fa):
                            head = sep.join(title.split(sep)[:-1])[1:]
                            count = title.split(sep)[-1]
                            header = ">{samp}.{id}_x{nr}\n".format(samp=sample, id=head, nr=count)
                            s = "{}\n".format(seq)
                            merged.write(header)
                            merged.write(s)
                            uniq_reads += 1
                            all_reads += int(count)
                            #zliczanie dlugosci wszystkich odczytów
                            for _ in range(int(count)):
                                lengths.append(int(len(seq.strip())))
                            #zliczanie dlugosci unikalnych odczytów
                            #lengths.append(int(len(seq.strip())))
                    lengths = np.array(lengths)
                    reads_count[in_file] = [uniq_reads, all_reads, np.round(np.mean(lengths),0)]
                else:
                    path, name = os.path.split(in_file)
                    if encoding == 'gzip':
                        coll_name = ".".join(name.split(".")[:-2]) + '.collapsed.fasta'
                    else:
                        coll_name = ".".join(name.split(".")[:-1]) + '.collapsed.fasta'
                    coll_path = os.path.join(path, coll_name)
                    coll_dict = self.collapse(in_file, coll_path, form, save, _open)
                    for s in coll_dict:
                        header = ">{samp}.{id}_x{nr}\n".format(samp=sample,id=coll_dict[s][1],nr=coll_dict[s][0])
                        seq = "{seq}\n".format(seq=s)
                        merged.write(header)
                        merged.write(seq)
                        uniq_reads += 1
                        all_reads += int(coll_dict[s][0])
                        #zliczanie dlugosci wszystkich odczytów
                        for _ in range(int(coll_dict[s][0])):
                            lengths.append(int(len(seq.strip())))
                        #zliczanie dlugosci unikalnych odczytów
                        #lengths.append(int(len(seq.strip())))
                    lengths = np.array(lengths)
                    reads_count[in_file] = [uniq_reads, all_reads, np.round(np.mean(lengths),0)]
        merged.close()
        formats = PrettyTable()
        formats.field_names = ["Format", "Number of files"]
        for form in form_dict:
            formats.add_row([form, form_dict[form]])
        formats.sortby = "Number of files"
        formats.align["Format"] = "l"
        formats.align["Number of files"] = "r"
        
        if plainTab:
            formats.set_style(PLAIN_COLUMNS)
        else: 
            formats.title = "Input files formats statistics"
        reads = PrettyTable()
        reads.field_names = ["File name", "Nr of unique reads", "Nr of all reads", "Mean length of reads"]
        for name in reads_count:
            c = reads_count[name]
            reads.add_row([name, c[0], c[1], c[2]])
        reads.sortby = "File name"
        reads.align["File name"] = "l"
        reads.align["Nr of uniq reads"] = "r"
        reads.align["Nr of all reads"] = "r"
        reads.align["Mean length of reads"] = "r"
        
        if plainTab:
            reads.set_style(PLAIN_COLUMNS)
            to_print = "{0}\n{1}\n{2}\n{3}\n".format("Input files formats statistics", formats.get_string(), "Input files reads statistics", reads.get_string())
        else:
            reads.title ="Input files reads statistics"
            to_print = "{0}\n{1}\n".format(formats.get_string(), reads.get_string())
        #print("Merged file name: ", out_file)
        print(to_print)
        if save_stat:# TODO -> do pliku - zawsze czy na życzenie użytkownika ? wypluć na ekran zawsze
            if gz:
                name = ".".join(out_file.split(".")[:-2]) + "_statistics.txt"
            else:
                name = ".".join(out_file.split(".")[:-1]) + "_statistics.txt"
            stat_file = open(name, "w")
            stat_file.write(to_print)
            stat_file.close()

    def exception_handling(self, msg):
        """When exception happend, print messeage and exit program

        Arguments:
            msg {string} -- Messeage to print
        """
        sys.stderr.write(msg)
        sys.exit()

    def collapse(self, inp: str, out: str, ext: str, save: bool, _open):
        """Collapse reads if it's needed

        Arguments:
            inp {str} -- input file with reads
            out {str} -- name of output file with collapsed reads
            ext {str} -- extension of file / file format: fasta/fastq
            save {bool} -- if true save collapsed file
            _open {[type]} -- open function for input file, regular open or gzip.open

        Returns:
            fa_dict -- dictionary with collapsed reads from input file
        """
        #cmd = 'fastx_collapser -i {input} -o {output}'.format(input=inp, output=out)
        #os.system(cmd)
        fa_dict = {}
        with _open(inp) as input_file:
            c=1 # ID sekwencji
            if ext == "fasta":
                for title, seq in SimpleFastaParser(input_file):
                    if seq not in fa_dict:
                        fa_dict[seq] = [1,c] # ile odczytów i jakie ID
                        c+=1
                    else:
                        fa_dict[seq][0] += 1
            elif ext == "fastq":
                for title, seq, qual in FastqGeneralIterator(input_file):
                    if seq not in fa_dict:
                        fa_dict[seq] = [1,c]
                        c+=1
                    else:
                        fa_dict[seq][0] += 1
        if save:
            coll = open(out,"w")
            for seq in fa_dict:
                header = ">{id}_x{count}\n".format(id=fa_dict[seq][1],count=fa_dict[seq][0])
                coll.write(header)
                coll.write(seq+"\n")
        return fa_dict

    def is_file(self, path: str):
        """Check if file exists

        Arguments:
            path {str} -- path to file

        Returns:
            [bool] -- True if file exists, if not print messeage and exit program
        """
        if not os.path.isfile(path):
            msg = "File error: No such file or directory -> {0}\n".format(path)
            self.exception_handling(msg)
        else:
            return True

    def is_empty(self, path: str):
        """Check if file is empty

        Arguments:
            path {str} -- path to file

        Returns:
            [bool] -- False if file has content, if not print messeage and exit program
        """
        if os.path.getsize(path) <= 0:
            msg = "File error: {0} -> file is empty\n".format(path)
            self.exception_handling(msg)
        else:
            return False
