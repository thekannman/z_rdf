//Copyright (c) 2015 Zachary Kann
//
//Permission is hereby granted, free of charge, to any person obtaining a copy
//of this software and associated documentation files (the "Software"), to deal
//in the Software without restriction, including without limitation the rights
//to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//copies of the Software, and to permit persons to whom the Software is
//furnished to do so, subject to the following conditions:
//
//The above copyright notice and this permission notice shall be included in all
//copies or substantial portions of the Software.
//
//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//SOFTWARE.

// ---
// Author: Zachary Kann

#include "xdrfile_trr.h"
#include "boost/program_options.hpp"
#include "z_sim_params.hpp"
#include "z_vec.hpp"
#include "z_constants.hpp"
#include "z_conversions.hpp"
#include "z_molecule.hpp"
#include "z_subsystem_group.hpp"
#include "z_histogram.hpp"
#include "z_gromacs.hpp"

namespace po = boost::program_options;
// Units are nm, ps.

int main (int argc, char *argv[]) {
  int st;
  SimParams params;

  enum RdfWeighting { kNone, kTotalDipole, kPermanentDipole, kInducedDipole };

  po::options_description desc("Options");
  desc.add_options()
    ("help,h",  "Print help messages")
    ("group1", po::value<std::string>()->required(),
     "Group around which rdf is centered")
    ("group2", po::value<std::string>()->required(),
     "Group for which rdf will be calculated around the positions of group 1")
    ("index,n", po::value<std::string>()->default_value("index.ndx"),
     ".ndx file containing atomic indices for groups")
    ("gro", po::value<std::string>()->default_value("conf.gro"),
     ".gro file containing list of atoms/molecules")
    ("top", po::value<std::string>()->default_value("topol.top"),
     ".top file containing atomic/molecular properties")
    ("output,o", po::value<std::string>()->default_value("rdf.txt"),
     "Name for output file.")
    ("max_time,t",
     po::value<double>()->default_value(0.0),
     "Maximum simulation time to use in calculations")
    ("rdf_weighting,w", po::value<std::string>()->default_value("none"),
     "weighting factor for rdf. Choose from 'none', 'dipole', "
     "'permanent_dipole', and 'induced_dipole'");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);
  if (vm.count("help")) {
    std::cout << desc << "\n";
    exit(EXIT_SUCCESS);
  }

  RdfWeighting rdf_weighting;
  bool need_field = false;
  if (vm["rdf_weighting"].as<std::string>() == "none") {
    rdf_weighting = kNone;
  } else if (vm["rdf_weighting"].as<std::string>() == "dipole") {
    rdf_weighting = kTotalDipole;
    need_field = true;
  } else if (vm["rdf_weighting"].as<std::string>() == "permanent_dipole") {
    rdf_weighting = kPermanentDipole;
  } else if (vm["rdf_weighting"].as<std::string>() == "induced_dipole") {
    rdf_weighting = kInducedDipole;
    need_field = true;
  }

  std::map<std::string, std::vector<int> > groups;
  groups = ReadNdx(vm["index"].as<std::string>());
  std::vector<Molecule> molecules = GenMolecules(vm["top"].as<std::string>(),
                                                 params);

  SystemGroup all_atoms(vm["gro"].as<std::string>(), molecules);
  SubsystemGroup *first_group_pointer =
      SubsystemGroup::MakeSubsystemGroup(
          vm["group1"].as<std::string>(),
          SelectGroup(groups, vm["group1"].as<std::string>()), all_atoms);
  SubsystemGroup &first_group = *first_group_pointer;
  SubsystemGroup *second_group_pointer =
      SubsystemGroup::MakeSubsystemGroup(
          vm["group2"].as<std::string>(),
          SelectGroup(groups, vm["group2"].as<std::string>()), all_atoms);
  SubsystemGroup &second_group = *second_group_pointer;

  rvec *x_in = NULL;
  matrix box_mat;
  arma::rowvec box = arma::zeros<arma::rowvec>(DIMS);
  std::string xtc_filename = "prod.xtc";
  XDRFILE *xtc_file;
  params.ExtractTrajMetadata(strdup(xtc_filename.c_str()), (&x_in), box);
  xtc_file = xdrfile_open(strdup(xtc_filename.c_str()), "r");
  params.set_box(box);
  params.set_max_time(vm["max_time"].as<double>());

  if (need_field)
    second_group.OpenFieldFile();

  Histogram rdf_hist(params.box(0)/2.0/0.05, 0.05);

  arma::rowvec dx, dipole;
  float time, prec;
  int step = 0;
  double avg_volume;
  for (step = 0; step < params.max_steps(); ++step) {
    if(read_xtc(xtc_file, params.num_atoms(), &st, &time, box_mat, x_in, &prec))
      break;
    params.set_box(box_mat);
    avg_volume += params.volume();

    first_group.set_positions(x_in);
    second_group.set_positions(x_in);

    if (need_field) {
      second_group.SetElectricField(all_atoms, params.box());
      if (!second_group.field_check())
        second_group.WriteElectricField();
    }
    second_group.UpdateCom();
    for (int i_other = 0; i_other < second_group.num_molecules(); ++i_other) {
      switch(rdf_weighting) {
        case kPermanentDipole:
          second_group.PermanentDipole(i_other, params.box(), dipole, true);
          dipole *= ENM_TO_D;
          break;
        case kInducedDipole:
          second_group.InducedDipole(i_other, dipole, true);
          dipole *= ENM_TO_D;
          break;
        case kTotalDipole:
          second_group.PermanentDipole(i_other, params.box(), dipole, true);
          second_group.InducedDipole(i_other, dipole);
          dipole *= ENM_TO_D;
          break;
        default:
          break;
      }
      for (int i_atom = 0; i_atom < first_group.size(); ++i_atom) {
        FindDxNoShift(dx, first_group.position(i_atom),
                      second_group.com_position(i_other), params.box());
        double distance = arma::norm(dx);
        double weight;
        if (rdf_weighting == kTotalDipole || rdf_weighting == kInducedDipole ||
            rdf_weighting == kPermanentDipole) {
          dx = arma::normalise(dx);
          weight = arma::dot(dx,dipole);
        } else {
          weight = 1.0;
        }
        rdf_hist.Add(distance, weight);
      }
    }

  }
  xdrfile_close(xtc_file);

  avg_volume /= step;
  rdf_hist.Multiply(1.0/avg_volume/step/first_group.size()/
                    second_group.num_molecules());
  rdf_hist.WeightByShellVolume();

  rdf_hist.Print(vm["output"].as<std::string>(), true);
}
 // main
