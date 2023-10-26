#ifdef HAVE_CONFIG_H
# include "config.h"     
#endif
#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>
#include <limits>

#include "dune/common/parallel/mpihelper.hh"
#include <dune/common/exceptions.hh> 

#include <dune/grid/uggrid.hh>  
#include <dune/grid/common/gridinfo.hh> 
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/common/fvector.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

// Izračunaj kut između p1-p2 i p3-p2
template <typename Point>
double kut(Point const & p1, Point const & p2, Point const & p3)
{
    Dune::FieldVector<double, 2> vector1 = p1 - p2;
    Dune::FieldVector<double, 2> vector2 = p3 - p2;

    // kut = cos^-1(a*b / (norm(a) * norm(b)))

    // Skalarni produkt
    double dotProd = dot(vector1,vector2);
    // Norme
    double normProd = vector1.two_norm() * vector2.two_norm();

    // Racunanje kuta
    double kut = std::acos(dotProd / normProd);

    // Prebacivanje u stupnjeve
    kut = 180 * kut / M_PI;

    return kut;
}

int main(int argc, char** argv)
{
    const int dim = 2;
    using GridType = Dune::UGGrid<dim>;
    using LeafGridView = GridType::LeafGridView;

    Dune::MPIHelper::instance(argc, argv);

    if(argc < 3){
       std::cerr << "Usage: " << argv[0]
                 << " ime_grid_datoteke.msh br_profinjenja" << std::endl;
       std::exit(1);
    }
    
    // gridView JE LeafGridView!
    double min_kut = 1000;
    double max_kut = -1000;
    int min_ind = 0;
    int max_ind = 0;
    int i = 0;
    double kut1;



    bool verbosity = true;
    bool insertBoundarySegments = false;  // Bez toga Dune::GmshReader zna podbaciti (u 3D)

    std::unique_ptr<GridType> pgrid = Dune::GmshReader<GridType>::read(argv[1],
                                                        verbosity, insertBoundarySegments);

    int no_r = std::stoi(argv[2]);  // broj profinjenja
    pgrid->globalRefine(no_r);     // profini mrežu
    auto gv = pgrid->leafGridView();
    for(auto const & element : elements(gv))
    {
        // Mreza se sastoji od trokuta, tj element.geometry().corners() == 3
        kut1 = kut(element.geometry().corner(0),element.geometry().corner(1), element.geometry().corner(2));
        if(kut1 < min_kut)
        {
            min_kut = kut1;
            min_ind = i;
        }
        if(kut1 > max_kut)
        {
            max_kut = kut1;
            max_ind = i;
        }

        kut1 = kut(element.geometry().corner(1),element.geometry().corner(2), element.geometry().corner(0));
        if(kut1 < min_kut)
        {
            min_kut = kut1;
            min_ind = i;
        }
        if(kut1 > max_kut)
        {
            max_kut = kut1;
            max_ind = i;
        }

        kut1 = kut(element.geometry().corner(2),element.geometry().corner(0), element.geometry().corner(1));
        if(kut1 < min_kut)
        {
            min_kut = kut1;
            min_ind = i;
        }
        if(kut1 > max_kut)
        {
            max_kut = kut1;
            max_ind = i;
        }
        i++;
    }


   // ISPISATI BROJ ELEMENATA; MINIMALNI I MAKSIMALNI KUT U STUPNJEVIMA:
    std::cout << "Broj elemenata: " << i << "\n";
    std::cout << "Minimalni kut: " << min_kut << " indeks elementa na kojem se postize: " << min_ind << "\n";
    std::cout << "Maksimalni kut: " << max_kut << " indeks elementa na kojem se postize: " << max_ind << "\n";

    // Ispis mreže u VTK formatu (u datoteku poluvijenac.vtu)
    Dune::VTKWriter<LeafGridView> vtkwriter(gv);
    vtkwriter.write("poluvijenac");

    return 0;
}
