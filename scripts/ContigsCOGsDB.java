
import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
public class ContigsCOGsDB {

	@SuppressWarnings("deprecation")
	public static void main(String[] args) throws IOException {
				FileReader rpsblastFile = null;
				FileWriter COGsFile = null;
				String xmlLine = "";
				
				URL u;
			    InputStream is = null;
			    DataInputStream dis;
			    String dir = System.getProperty("user.dir");
			    
			    long startTime = System.currentTimeMillis();
				try {
					try {
						//COGsFile = new FileWriter(dir+"/Output/PROKKA_COGsDB.tsv"); //normally results should go to /Output folder
						COGsFile = new FileWriter(dir+"/PROKKA_COGsDB.tsv"); //due to memory shortage though will store them to hard drive
						COGsFile.append("CDS_ID\tCOG\tName\tStart\tEnd\n");
						
						rpsblastFile = new FileReader(new File(args[0]));
						@SuppressWarnings("resource")
						BufferedReader rpsblastSc = new BufferedReader(rpsblastFile);

						int cdsNum = 0;
						String line;
						while((line = rpsblastSc.readLine()) != null) {
							String[] lineTok = line.split("\\t", -1);
							String idCDD = lineTok[1].split("\\|")[2];
							
							try{
								u = new URL("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=cdd&id="+idCDD);
								is = u.openStream();
								dis = new DataInputStream(new BufferedInputStream(is));
	
								String COG = "", subtitle = "";
								while ((xmlLine = dis.readLine()) != null) {
									if(xmlLine.contains("Accession")) {
										COG = xmlLine.split(">")[1].replace("</Item", "");
										 //COG = StringUtils.substringBetween(xmlLine, ">", "<");
									} else if(xmlLine.contains("Subtitle")) {
										subtitle = xmlLine.split(">")[1].replace("</Item", "");
										//subtitle = StringUtils.substringBetween(xmlLine, ">", "<");
										break;
									}
								}		

								if((cdsNum % 10000)==0) {
									System.out.println(cdsNum);
								}
								cdsNum ++;
						        COGsFile.append(lineTok[0]+"\t"+COG+"\t"+subtitle+"\t"+lineTok[5]+"\t"+lineTok[6]+"\n");
						        
						        dis.close();
						        is.close();
					        } catch (IOException ex) {
					        	//System.err.println("File exception");
					        }
						}
					} finally {
						rpsblastFile.close();
						COGsFile.close();
					}
		        } catch (IOException ex) {
		        	System.err.println("File exception");
		        }
				
				long endTime = System.currentTimeMillis();
		   		System.out.println("Excution took " + (endTime - startTime) + " milliseconds");
	}

}
