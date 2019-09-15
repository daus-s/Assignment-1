#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <string.h>
#include <stdio.h>
#include <cstdio>
#include <stdlib.h>
#include <ctype.h>
#include <time.h>

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

        double length(double x, double y);
        string generate(double x, double y);
        void parse(string data);
        double mean(string data);
        double stddev(string data);
        int sum();
        double variance(string data);
        double gaussian(double x, double y);

        DNA();
        DNA(string data);
        ~DNA();
};

int main(int argc, char** argv)
{
    cout << "dna analyzer - dausMcChar" << endl;
    //checks to see if there is an arguement passed into main from terminal
    if (argc > 1)
    {
        //variable to check consent of user *******************
        char cont = 'y';
        while (true)
        {
            //pathname is name of file inputted
            char* pathname = argv[1];
            cout << "opening "<< pathname << "..." << endl;
            ifstream analyze;
            analyze.open(pathname, ifstream::in);
            // if (analyze.is_open())
            // {
            //     cout << "open" << endl;
            // }
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
            while (analyze.good())
            {
                if (line.back() == '\n')
                    line = line.substr(0, line.size()-2); //ejfibsaojdnflksadnflkdsmnalkfnsdkljnfkmsngkm sdkmgnksdjngkmdsn gkn sdkg sdkn gvkn there is probably an error here

                file += line + " ";
                getline(analyze, line);


                //cout << "line: " << line << endl;
            }
            cout << "finished reading, added to string file" << endl;
            analyze.close();
            DNA dna = DNA(file);

            //F in the chat boys, closes ifstream because no need at this time


            //Analysis time babey what it do
            //debug lines
            cout << "appending to daus.out"  << endl;

            clout << "probability of nucleotides" << endl;
            clout << "\ta:" << ((double)dna.a/dna.nucleotides) << endl;
            clout << "\tc:" << ((double)dna.c/dna.nucleotides) << endl;
            clout << "\tg:" << ((double)dna.g/dna.nucleotides) << endl;
            clout << "\tt:" << ((double)dna.t/dna.nucleotides) << endl;

            clout << "propability of bigrams:" << endl;
            clout << "\taa:" << ((double)dna.aa/(dna.nucleotides-1)) << endl;
            clout << "\tac:" << ((double)dna.ac/(dna.nucleotides-1)) << endl;
            clout << "\tag:" << ((double)dna.ag/(dna.nucleotides-1)) << endl;
            clout << "\tat:" << ((double)dna.at/(dna.nucleotides-1)) << endl;
            clout << "\tca:" << ((double)dna.ca/(dna.nucleotides-1)) << endl;
            clout << "\tca:" << ((double)dna.cc/(dna.nucleotides-1)) << endl;
            clout << "\tca:" << ((double)dna.cg/(dna.nucleotides-1)) << endl;
            clout << "\tct:" << ((double)dna.ct/(dna.nucleotides-1)) << endl;
            clout << "\tga:" << ((double)dna.ga/(dna.nucleotides-1)) << endl;
            clout << "\tgc:" << ((double)dna.gc/(dna.nucleotides-1)) << endl;
            clout << "\tgg:" << ((double)dna.gg/(dna.nucleotides-1)) << endl;
            clout << "\tgt:" << ((double)dna.gt/(dna.nucleotides-1)) << endl;
            clout << "\tta:" << ((double)dna.ta/(dna.nucleotides-1)) << endl;
            clout << "\ttc:" << ((double)dna.tc/(dna.nucleotides-1)) << endl;
            clout << "\ttg:" << ((double)dna.tg/(dna.nucleotides-1)) << endl;
            clout << "\ttt:" << ((double)dna.tt/(dna.nucleotides-1)) << endl;

            clout << "standard deviation: " << dna.stddev(file) << endl;
            clout << "variance: " << dna.variance(file) << endl;
            clout << "mean: " << dna.mean(file) << endl;
            clout << "sum: " << dna.sum() << endl << endl;


            for (int i = 0; i < 1000; i++)
            {
                double x = (double)rand() / RAND_MAX;
                double y = (double)rand() / RAND_MAX;
                clout << dna.generate(x, y) << endl;
            }
            //checks next char to see if user wants to keep going
            cout << "Continue? (y/n)" << endl;
            cin >> cont;
            if (cont == 'n' || cont == 'N')
            {
              return 0;
              //closes off output stream to daus.out
              clout.close();
            }
        }
    }
    return 2;
}

//returns stanardly distributed length
double DNA::gaussian(double x, double y)
{
    double z = sqrt(-2*log(x)) * cos(2*M_PI*y);
    return abs(z);
}

double DNA::length(double x, double y)
{
    //cout << "gaussian: " << gaussian(x,y) << "mean: " << mean(sequences) << endl;
    return (stddev(sequences) * gaussian(x,y)) + mean(sequences);
}

string DNA::generate(double x, double y)
{
    int size = (int)length(x,y);
    //cout << size << endl;
    string deoxy = "";
    int numA = size * ((double)a/nucleotides);
    int numC = size * ((double)c/nucleotides);
    int numG = size * ((double)g/nucleotides);
    int numT = size * ((double)t/nucleotides);
    //size = numA + numC + numG + numT;

    for (int i = 0; i < size; i++)
    {
        int n = rand() % 4;
        if (n == 0)
        {
            if (numA != 0)
            {
              deoxy += 'a';
              numA--;
            }
            else
            {
                n++;
            }
        }
        if (n == 1)
        {
            if (numC != 0)
            {
              deoxy += 'c';
              numC--;
            }
            else
            {
                n++;
            }
        }
        if (n == 2)
        {
            if (numG != 0)
            {
              deoxy += 'g';
              numG--;
            }
            else
            {
                n++;
            }
        }
        if (n == 3)
        {
            if (numT != 0)
            {
              deoxy += 't';
              numT--;
            }
            else
            {
                n = 0;
            }
        }
    }

    return deoxy;

}

//this is called A$AP when 'file' is initialized
void DNA::parse(string data)
{

      string index = " ";
      for (int i = 0; i < data.size(); ++i)
      {
          index = data[i];
          if (index == "a")
              a++;
          if (index == "c")
              c++;
          if (index == "g")
              g++;
          if (index == "t")
              t++;
          if (index == " ")
              strands++;

          if (index == " "|| data[i+1] == ' ')
          {
             continue;
          }
          //cout << "bigram " << i << ":" << index << data[i+1] << endl; //this one works the other one (below) doesnt
          string bigram = index + data[i+1];
          if (i < data.size()-1 && bigram.compare("aa")  == 0)
              aa++;
          if (i < data.size()-1 && bigram.compare("ac")  == 0)
              ac++;
          if (i < data.size()-1 && bigram.compare("ag")  == 0)
              ag++;
          if (i < data.size()-1 && bigram.compare("at")  == 0)
              at++;
          if (i < data.size()-1 && bigram.compare("ca")  == 0)
              ca++;
          if (i < data.size()-1 && bigram.compare("cc")  == 0)
              cc++;
          if (i < data.size()-1 && bigram.compare("cg")  == 0)
              cg++;
          if (i < data.size()-1 && bigram.compare("ct")  == 0)
              ct++;
          if (i < data.size()-1 && bigram.compare("ga")  == 0)
              ga++;
          if (i < data.size()-1 && bigram.compare("gc")  == 0)
              gc++;
          if (i < data.size()-1 && bigram.compare("gg")  == 0)
              gg++;
          if (i < data.size()-1 && bigram.compare("gt")  == 0)
              gt++;
          if (i < data.size()-1 && bigram.compare("ta")  == 0)
              ta++;
          if (i < data.size()-1 && bigram.compare("tc")  == 0)
              tc++;
          if (i < data.size()-1 && bigram.compare("tg")  == 0)
              tg++;
          if (i < data.size()-1 && bigram.compare("tt")  == 0)
              tt++;

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
    for (int i = 0; i < data.size(); ++i)
    {
        tolower(data[i]);
    }
    parse(data);
    sequences = data;
}
DNA::DNA()
{}
DNA::~DNA()
{}
