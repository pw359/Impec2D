#ifndef GNUPLOTTER_H_
#define GNUPLOTTER_H_

#include "SimParams.h"
#include "GenSimDefs.h"
class Gnuplot;

class GnuPlotter {

	Gnuplot * p;
	string cmds;
	bool firstplot;
	string temp_dir;

public:
	GnuPlotter(bool levelset, const string temp_dir);
	GnuPlotter(const string temp_dir);
	virtual ~GnuPlotter();
	void Plot(SimParams const * const sp, const Array2d& arr);
	void ContourPlot(const string & data);
	void PlotConvergence(const string& data);
	void ReplaceTokenByFile(string & s, const string & token, const string & file);
};

#endif /* GNUPLOTTER_H_ */
