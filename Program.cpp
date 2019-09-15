#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <string.h>
#include <stdio.h>


using namespace std;

class DNA
{
    //count of each bigram and nucleotide
    public:
        int aa = 0;
        int ac = 0;
        int ag = 0;
        int at = 0;
        int ca = 0;
        int cc = 0;
        int cg = 0;
        int ct = 0;
        int ga = 0;
        int gc = 0;
        int gg = 0;
        int gt = 0;
        int ta = 0;
        int tc = 0;
        int tg = 0;
        int tt = 0;
        //nucleotides
        int a = 0;
        int c = 0;
        int g = 0;
        int t = 0;
        int nucleotides = 0; //this is everysingle nucleotide, ==a+c+g+t
        int strands = 0; //this is the count of lines
        string sequences = "";

        int length(double a, double b);
        string generate();
        void parse(string data);
        double mean(string data);
        double stddev(string data);
        int sum();

        DNA();
        DNA(string data);
        ~DNA();
};

int main(int argc, char** argv)
{
    cout << "dna" << endl;
    //checks to see if there is an arguement passed into main from terminal
    if (argc > 1)
    {
        //variable to check consent of user *******************
        char cont = 'y';
        while (cont != 'n' || cont != 'N')
        {
            //pathname is name of file inputted
            char* pathname = argv[1];
            ifstream analyze;
            analyze.open(pathname, ifstream::in);
            string line;

            //writes line from file to new varaible representing a line
            getline(analyze,line);
            //readline is java, getline is c++

            ofstream clout;
            //compiled line output, also https://twitter.com/dausMcChar/status/1168976276046143488?s=20
            clout.open("daus.out", ofstream::trunc);
            string file = "";
            clout << "Daus(Davis) Carmichael \n2328047\ncarmichael@chapman.edu" << endl;

            //this loop reads file and inputs into my file variable
            //checks to see if the new line is end of file, if it is we stop analysis
            while (!analyze.eof())
            {
                if (line.back() == '\r')
                    line = line.substr(0, line.size()-2); //ejfibsaojdnflksadnflkdsmnalkfnsdkljnfkmsngkm sdkmgnksdjngkmdsn gkn sdkg sdkn gvkn there is probably an error here
                file += line + " ";
                getline(analyze, line);
            }
            analyze.close();
            DNA dna = DNA(file);

            //F in the chat boys, closes ifstream because no need at this time


            //Analysis time babey what it do




            //checks next char to see if user wants to keep going
            cout << "Continue? (y/n)" << endl;
            cin >> cont;
        }
    }
    return 2;
}

//returns stanardly distributed length
int DNA::length(double x, double y)
{
    int z = sqrt(-2*log(x)) * cos(2*M_PI*y);
    return abs(z);
}

string DNA::generate()
{

}

//this is called A$AP when 'file' is initialized
void DNA::parse(string data)
{
      char index = ' ';
      for (int i = 0; i < data.size(); ++i)
      {
          index = data[i];
          if (index == 'a')
              a++;
          if (index == 'c')
              c++;
          if (index == 'g')
              g++;
          if (index == 't')
              t++;
          if (index == ' ')
              strands++;
      }
    nucleotides = a + c + g + t;
}
//returns the mean number of nucleotides per dna strand
double DNA::mean(string data)
{
    return (double)nucleotides / (double)strands;
}

int DNA::sum()
{
    return nucleotides;
}

double DNA::stddev(string data)
{
    int strandLength;
    double sum = 0.0;
    double avg = mean(data);

    for (int i = 0; i < nucleotides; ++i)
    {
        if (data[i] != ' ')
        {
            strandLength++;
        }
        if (data[i]==' ')
        {
            sum += pow(strandLength - avg, 2);
        }
    }

    return sqrt(sum/strands);
}
double DNA::variance(string data)
{
  return pow(stddev(data),2);
}
//constructs DNA object
DNA::DNA(string data)
{
    parse(data);
    sequences = data;
}
DNA::DNA()
{}
DNA::~DNA()
{}
