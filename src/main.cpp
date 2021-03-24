/// @brief: This program is a small utility which reads a list of numbers
///         from one or more text files and prints out the correlation function.


#include <vector>
#include <cassert>
#define _USE_MATH_DEFINES
#include <cmath>
using namespace std;
#include "io.hpp"
#include "ndautocrr.hpp"
#include "err.hpp"



// global variables
const char g_program_name[]   = "ndautocrr";
const char g_version_string[] = "0.12.1";
const char g_date_string[]    = "<2021-3-24>";
const char *g_spaces = " \t"; // whitespace characters (exluding newlines)
const char *g_spaces_and_newlines = " \t\n";





int
main(int argc, char **argv)
{
  try {
    cout.precision(14);

    cerr
      << g_program_name   << ", v"
      << g_version_string << " "
      << g_date_string    << "\n";


    bool is_periodic = false;
    size_t L = 0; // an impossible starting value

    // Now parse the argument list:

    double threshold = -1.01; //(setting this to a value below -1 disables)
    bool subtract_ave = true;
    bool report_rms = false;
    bool report_nsum = false;

    {
      bool syntax_error_occured = false;

      long i=1;
      while (i < argc)
      {
        int ndelete = 0;
        if (strcmp(argv[i], "-L")==0) {
          if ((strlen(argv[i+1]) == 0) || (! isdigit(argv[i+1][0])))
            throw InputErr("Error: Expected a number following the -L flag.\n");
          else {
            L = atoi(argv[i+1]);
          }
          ndelete = 2;
        }
        else if ((strcmp(argv[i], "-p")==0) ||
                 (strcmp(argv[i], "-P")==0) ||
                 (strcmp(argv[i], "-periodic")==0) ||
                 (strcmp(argv[i], "--periodic")==0))
        {
          is_periodic = true;
          cerr <<"Argument found: "<<argv[i]<<" -> PERIODIC BOUNDARY CONDITIONS USED.\n" << endl;
          ndelete = 1;
        }
        else if (strcmp(argv[i], "-ave")==0)
        {
          subtract_ave = true;
          ndelete = 1;
        }
        else if (strcmp(argv[i], "-avezero")==0)
        {
          subtract_ave = false;
          ndelete = 1;
        }
        else if (strcmp(argv[i], "-rms")==0)
        {
          report_rms = true;
          ndelete = 1;
        }
        else if (strcmp(argv[i], "-nsum")==0)
        {
          report_nsum = true;
          ndelete = 1;
        }
        else if ((strcmp(argv[i], "-t")==0) ||
                 (strcmp(argv[i], "-T")==0) ||
                 (strcmp(argv[i], "-threshold")==0) ||
                 (strcmp(argv[i], "--threshold")==0))
        {
          syntax_error_occured = false;

          if ((argc > i+1) &&
              (isdigit(argv[i+1][0]) ||
               (argv[i+1][0] == '-') ||
               (argv[i+1][0] == '+') ||
               (argv[i+1][0] == '.') ||
               (argv[i+1][0] == 'e') ||
               (argv[i+1][0] == 'E')))
            threshold = atof(argv[i+1]);
          else
            syntax_error_occured = true;

          if ((threshold < -1.0) || (threshold > 1.0))
            syntax_error_occured = true;

          if (syntax_error_occured)
          {
            cerr <<
              "Error: Expected a number between -1.0 and 1.0 following the -t flag.\n";
            if (argc > i+1)
            {
              stringstream err_msg;
              err_msg <<
                "       (This \"threshold\" should be expressed as a fraction of <(x-<x>)^2>)\n"
                "       Instead, you specified \""<<argv[i]<<
                " "<<argv[i+1]<<"\"\n" << endl;
              throw InputErr(err_msg.str().c_str());
            }
          }

          cerr <<
            "The correlation function will stop when dropping below a threshold.\n"
            "threshold = " << threshold << " (relative to the peak at separation 0).\n";
          ndelete = 2;
        }

        if (ndelete > 0) { // if the argument was recogized
          // Delete the argument(s) recognized in this pass
          for (long j = i; j+ndelete < argc; j++)
            argv[j] = argv[j+ndelete];
          argc -= ndelete;
        }
        else
          i++; //we will deal with unrecongized arguments later

      } // while (i < argc)
    } // command line argument parsing ends here


    // Any left over arguments?  Hope not.
    if (argc > 1) {
      stringstream err_msg;
      err_msg << "Unnexpected argument: \"" << argv[1] << "\"\n";
      throw InputErr(err_msg.str().c_str());
    }


    // allocate the array to store the auto-correlation function

    NdAutocrr<double>
      ndautocrr = NdAutocrr<double>(threshold,
                                    L,
                                    is_periodic,
                                    subtract_ave,
                                    report_rms);


    // now read in the data from the file

    long n_data_sets = 1;
    vector<double> vX_d;
    vector<vector<double> > vvX_id;
    g_filename.assign("standard-input/terminal");
    g_line=1;           //keep track of which line number
    Skip(cin, g_spaces_and_newlines);
    long prev_line = g_line; //used to figure out if 2 numbers on same line
    long L_min = -1;
    while(cin)
    {
      double x;
      //cin >> x;
      x = ReadScalar<double>(cin, g_spaces_and_newlines);
      if (! cin) break;
      assert(g_line == prev_line);
      vX_d.push_back(x);

      if (cin) Skip(cin, g_spaces_and_newlines);
      //did the line number increment or not?
      if (g_line > prev_line) {
        assert(vX_d.size() > 0);
        vvX_id.push_back(vX_d);
        vX_d.resize(0);
      }
      if (g_line - prev_line > 1) {
        if (vvX_id.size() > 0) {
          cerr << "#  processing data set #" << n_data_sets << endl;
          if ((n_data_sets > 1) && (threshold > -1.0)) {
            throw InputErr("ERROR: Do not use -threshold when analyzing files containing multiple data\n"
                           "       sets separated by blank lines (sometimes also called \"trajectories\").\n"
                           "       Use the -L argument instead.\n");
          }

          ndautocrr.AccumulateSingle(vvX_id);
          n_data_sets++;
        }
        vvX_id.resize(0);
      }
      prev_line = g_line;
    }

    if (vvX_id.size() > 0) {
      cerr << "#  processing data set #" << n_data_sets << endl;
      if ((n_data_sets > 1) && (threshold > -1.0)) {
        throw InputErr("ERROR: Do not use -threshold when analyzing files containing multiple data\n"
                       "       sets separated by blank lines (sometimes also called \"trajectories\").\n"
                       "       Use the -L argument instead.\n");
      }

      ndautocrr.AccumulateSingle(vvX_id);
      n_data_sets++;
    }


    ndautocrr.Finalize();


    //Now print the corrlation function to the standard out
    //assert(L <= vCsum.size());

    cerr << "#----- delta  C(delta) -----\n" << endl;

    L = ndautocrr.size();

    for (size_t j=0; j <= L; ++j)
    {
      if (ndautocrr.vNumSamples[j] > 0) {
        cout << j
             << " " << ndautocrr.vC[j];
        if (report_rms)
          cout << " " << ndautocrr.vCrms[j];
        if (report_nsum)
          cout << " " << ndautocrr.vNumSamples[j];
        cout << "\n";
      }
    }

    // Now print back the corrlation length
    double correlation_length = ndautocrr.GuessCorrelationLength();

    cerr <<
      "\n"
      "#--------------------------------------\n"
      "# correlation length = " << correlation_length
      // << "\n"
      //"# g = 1 + 2*(correlation length) = "
      //   << 1.0 + (2.0*correlation_length)
         << endl;

  } // try
  catch (const std::exception& e) {
    cerr << "\n" << e.what() << endl;
    exit(1);
  }
}


