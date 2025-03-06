# include <iostream>
# include <string>
# include <fstream>
# include <vector>

using namespace std;

void RENAME_CONTIGS(string main_fasta_file, string new_single_line_fasta_file)
{
	string line, seq_description, seq;
	fstream main_contig_file, new_fasta_file;
	main_contig_file.open(main_fasta_file, ios::in);
	new_fasta_file.open(new_single_line_fasta_file, ios::app);
	
	string prev_seq = "";
    long long cnt = 0;
	while (getline (main_contig_file, line)) { 
  		
  		if(line[0] == '>'){
  			//previous seq writing
  			
  			if(prev_seq.empty() || prev_seq.length() < 1000){ // filter out contigs with length less than 1000 bp
  				// do nothing
  			}
  			else{
                cnt++;
  				new_fasta_file << ">contig_" + to_string(cnt) + "\n";
				new_fasta_file << prev_seq + "\n";
  			}
  			
  			// now this seq
  			seq_description = line;
  			prev_seq = "";
  		}
  		else{
  			prev_seq+=line;
  			
  		}
		
	}
	// handling last sequence in the fasta file
    if(prev_seq.length() >= 1000){ // filter out contigs with length less than 1000 bp
        cnt++;
        new_fasta_file << ">contig_" + to_string(cnt) + "\n";
        new_fasta_file << prev_seq + "\n";
  	}
	
	return;
}


int main(int argc, char *argv[])
{	
	
	string corresponding_contig_file = string(argv[1]);
    string out_prefix = string(argv[2]);

	string new_contig_file = out_prefix + "_renamed_contigs.fa";
	
	RENAME_CONTIGS(corresponding_contig_file, new_contig_file);
	
	return 0;

}