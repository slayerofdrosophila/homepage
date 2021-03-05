`
// Name: Joshua Liu
// COSI 12b, Spring 2021
// Programming Assignment 2
// Description: Parses a file with DNA information on it. Each line pair has a region name and nucleotides.
// This program determines the properties of the region of nucleotides and outputs them to a file.
// The file parsed and written to are specified by the user.


import java.io.*;    
import java.util.*;

public class DNA {
	
	static int minimumcodons = 5;
	static int cg_percentage = 30;
	static int nucleotidesnumber = 4;
	static int nucspercodon = 3;
	
	public static void main(String[] args) throws FileNotFoundException {
		
		Scanner console  = new Scanner(System.in);
		System.out.print("Input file name? ");
		String filename = console.next();
		
		System.out.print("Output file name? ");
		String outputfilename = console.next();
		
		console.close();
		
		int linenum = 0;
		char[] nucleotides = null;
		String regionname = "";
		String uppercaseme = "";
		
		Scanner filescan = new Scanner(new File(filename));
		PrintStream output = new PrintStream(new File(outputfilename));
		
		while (filescan.hasNextLine()) { // stops when document ends 
			String line = filescan.nextLine();
			linenum++;
			
			if (linenum % 2 == 0) { // if even line number, process this line
				uppercaseme = line.toUpperCase();
				nucleotides = toarray(uppercaseme);
				
			}
			else { // if odd line number, it must be the region name
				regionname = line;
			}
			
			if (linenum % 2 == 0){ // after every 2 lines, by which both vars should have been changed, print the info stored
				printinfo(nucleotides, regionname, output);
			}
		}
	}
	
	public static void printinfo(char[] nucleotides, String regionname, PrintStream output) { // prints the stuff. yea.
		
		output.println("Region Name: " + regionname);
		
		char[] noDashArray = removedashes(nucleotides);
		String nucStringNoDashes = new String(noDashArray);
		output.println("Nucleotides: " + nucStringNoDashes); // done differently so it's all in a string w/ no brackets and commas
		
		String nucCountsForPrint = Arrays.toString(count(nucleotides, 1));
		output.println("Nuc. Counts: " + nucCountsForPrint);
		
		float[] massForPrintArray = mass(count(nucleotides, 0));
		output.println("Total Mass%: " + Arrays.toString(massForPrintArray));
		
		String[] codonlistarray =  codonlist(removedashes(nucleotides));
		output.println("Codons List: " + Arrays.toString(codonlistarray));
		
		output.println("Is Protein?: " + isprotein(codonlistarray, massForPrintArray));
		
		output.println();
	}
	
	// REMEMBER, ACGT ORDER
	
	public static String isprotein(String[] codons, float[] masses) { // params: codons list, mass% list
		
		if (codons[0].equals("ATG") == false) {
			return "NO";
		}
			
		int stopcounter = 0;
		if (codons[codons.length - 1].equals("TAA") == false) {
			stopcounter++;
		}
		if (codons[codons.length - 1].equals("TAG") == false) {
			stopcounter++;
		}
		if (codons[codons.length - 1].equals("TGA") == false) {
			stopcounter++;
		}
		if (stopcounter >= 3) { // if it fails all 3 stop codon checks
			return "NO";
		}
				
		if (codons.length < minimumcodons) {
			return "NO";
		}
		if (masses[1] + masses[2] < cg_percentage) {
			return "NO";
		}
		
		return "YES";
	}
	
	public static String[] codonlist(char[] input) { //returns string[] with each 3-char codon in each index
		String[] output = new String[input.length / 3]; // btw, must have input array with no dashes
		
		char c1 = ' ';
		char c2 = ' ';
		char c3 = ' ';
		int counter = 0;
		int codoncounter = 0;
		
		for (int i = 0; i < input.length; i++) {
			if (counter == 0) {
				c1 = input[i];
			}
			if (counter == 1) {
				c2 = input[i];
			}
			if (counter == 2) {
				c3 = input[i];
			}
			counter++;
			
			if (counter >= nucspercodon) { // every 3 characters, reset counter 
								// and add the chars to the next slot in the output array
				counter = 0;
				output[codoncounter] = "" + c1 + c2 + c3;
				codoncounter++;
			}
		}
		return output;
	}
	
	
	public static float[] mass(int[] input) { // ACGT ORDER
		
		double[] masses = {135.128, 111.103, 151.128, 125.107, 100.000};
		for (int i = 0; i < 5; i++) {
			masses[i] = masses[i] * 100;
			masses[i] = Math.round(masses[i]);
			masses[i] = masses[i] / 100;
		}
		
		float a = (float) (input[0] * masses[0]);
		float c = (float) (input[1] * masses[1]);
		float g = (float) (input[2] * masses[2]);
		float t = (float) (input[3] * masses[3]);
		float junk = (float) (input[4] * masses[4]);
		
		float sum = a + c + g + t + junk;
		
		float[] output = new float[nucleotidesnumber + 1];
		
		output [0] = ((a / sum)* 100);
		output [1] = ((c / sum)* 100);
		output [2] = ((g / sum)* 100);
		output [3] = ((t / sum)* 100);
		output [4] = ((junk / sum)* 100);
		
		return output;
	}
	
	public static int[] count(char[] dna, int mode) { // ACGT ORDER
		
		int a = 0;
		int c = 0;
		int g= 0 ;
		int t = 0;
		int junk = 0;
		
		int[] nucslist = {a,c,g,t,junk};
		
		for (int i = 0; i < dna.length; i++) {
			 if (dna[i] == 'A') {
				 a++;
			 }
			 if (dna[i] == 'C') {
				 c++;
			 }
			 if (dna[i] == 'G') {
				 g++;
			 }
			 if (dna[i] == 'T') {
				 t++;
			 }
			 if (dna[i] == '-') {
				 junk++;
			 }
		}
		
		if (mode == 1) { // mode 1 is for printing
			int[] results = new int[nucleotidesnumber]; // 4 default +++++++++++++++++++++++++++++++++++++++++++++++++=============
			for (int n = 0; n < nucleotidesnumber; n++) {
				results[n] = nucslist[n];
			}
			return results;
		}
		else { // mode 0 is for internal program use, includes junk count
			int[] results = new int[nucleotidesnumber + 1]; // 5 default, +1 is for junk
			results[0] = a;
			results[1] = c;
			results[2] = g;
			results[3] = t;
			results[4] = junk;
			return results;
		}
		
	}
	
	public static char[] removedashes(char[] input) {
		
		int dashcount = 0;
		for (int i = 0; i < input.length; i++) {
			if (input[i] == '-') {
				dashcount++;
			}
		}
		char[] output = new char[input.length - dashcount];
		
		int counter = 0;
		for (int i = 0; i < input.length; i++) {
			if (input[i] != '-') {
				output[counter] = input[i];
				counter++;
			}
		}
		return output;
	}
	
	public static char[] toarray(String input) {
		char[] array = new char[input.length()];
		for (int i = 0; i < input.length(); i++) {
			array[i] = input.charAt(i);
		}
		return array;
	}
	
}
`
