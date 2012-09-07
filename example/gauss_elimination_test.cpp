#include "../include/gauss_elimination.h"

using namespace std;
using namespace jtf::algorithm;

//  2x2 + 2x1 + 2x1 = 2
//         x1 = 1
//  5x1 + x2 = 2;

int main()
{
  typedef double val_type;
  vector<val_type> node(4,0);
  boost::dynamic_bitset<> node_flag(4);

  jtf::algorithm::gauss_eliminator<val_type> ge(node, node_flag);

  {
    equation<val_type> eq;
    eq.add_expression(jtf::algorithm::make_expression(static_cast<size_t>(1),
                                                      static_cast<val_type>(2)));
    eq.add_expression(jtf::algorithm::make_expression(static_cast<size_t>(0),
                                                      static_cast<val_type>(2)));
    eq.add_expression(jtf::algorithm::make_expression(static_cast<size_t>(0),
                                                      static_cast<val_type>(2)));
    eq.value() = 2;
    ge.add_equation(eq);
  }
  {
    equation<val_type> eq;
    eq.add_expression(jtf::algorithm::make_expression(static_cast<size_t>(0),
                                                      static_cast<val_type>(1)));
    eq.value() = 1;
    ge.add_equation(eq);
  }
  {
    equation<val_type> eq;
    eq.add_expression(jtf::algorithm::make_expression(static_cast<size_t>(0),
                                                      static_cast<val_type>(5)));

    eq.add_expression(jtf::algorithm::make_expression(static_cast<size_t>(1),
                                                      static_cast<val_type>(1)));

    eq.value() = 2;
    ge.add_equation(eq);
  }
  for(size_t t = 0; t < node_flag.size(); ++t){
    if(node_flag[t] == true){
      cerr << "# node " << t << " = " << node[t] << endl;
    }else
      cerr << "# node " << t << " unknown." << endl;
  }
  return 0;
}
