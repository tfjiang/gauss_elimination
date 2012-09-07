#ifndef GAUSS_ELIMINATION_H
#define GAUSS_ELIMINATION_H

#include <map>
#include <vector>
#include <list>
#include <deque>
#include <iostream>
#include <boost/unordered_map.hpp>
#include <boost/static_assert.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/dynamic_bitset.hpp>
//#include <zjucad/matrix/matrix.h>

namespace jtf{
namespace algorithm{


template <typename T>
/**
 * @brief store an expression A * i, where A is the coeffieient
 *
 */
class expression{
public:
  expression():index(-1),coefficient(0){}
  expression(const size_t & index_, const T & coefficient_)
    : index(index_), coefficient(coefficient_){
    if(static_cast<int>(index_) < 0)
      std::cerr << "# [error] this index is negative." << std::endl;
  }
  size_t index;
  T coefficient;
  /**
   * @brief operator < is used to determin which expression ha
   *
   * @param other input other expression
   * @return bool
   */
  bool operator < (const expression<T> & b) const{
    return index < b.index;
  }

  /**
   * @brief To determin whether coefficent of this expression is zeros
   *
   * @return bool
   */
  bool is_zero() const{ return fabs(coefficient) < 1e-6;}
};

/**
 * @brief make expression with node_idx and value
 *
 * @param node_idx
 * @param value
 * @return expression<T>
 */
template <typename T>
expression<T> make_expression(const size_t & node_idx, const T & value)
{
  expression<T> temp(node_idx, value);
  return temp;
}

/**
 * @brief equation class
 *
 */
template <typename T>
class equation{
public:
  typedef typename std::list<expression<T> >::const_iterator eq_const_iterator;
  typedef typename std::list<expression<T> >::iterator eq_iterator;

  equation():value_(static_cast<T>(0)){}

  eq_const_iterator begin() const {return e_vec_.begin();}
  eq_iterator begin(){return e_vec_.begin();}

  eq_const_iterator end() const {return e_vec_.end();}
  eq_iterator end(){return e_vec_.end();}

  /**
   * @brief used to standardizate the equation, sort the expressions and merge
   *        similar items and normalization to make the first item's coefficient
   *        equal 1
   *
   * @return int
   */
  int standardization(){
    sort_equation();
    merge_similar();
    normalization();
    return 0;
  }

  /**
   * @brief update the equation with given node idx and value
   *
   * @param node_idx input node idx
   * @param node_value input node value
   * @return int if nothing change return 0, or return changes number
   */
  int update(const size_t & node_idx, const  T & node_value);
  /**
   * @brief sort equation accronding to the expressions
   *
   * @return int
   */
  int sort_equation(){
    e_vec_.sort();
    return 0;
  }

  /**
   * @brief merge similar items, WARNING. this function should be used after
   *        sorting.
   *
   * @return int
   */
  int merge_similar();

  /**
   * @brief normalize the equations: to make the first expression's coefficient
   *        equals 1
   *
   * @return int
   */
  int normalization();

  /**
   * @brief get the state of equation:
   *            if there are no expressions,
   *               if value != 0 return -1; error equation
   *               else return 0; cleared
   *            else
   *               if there is one expression return 1; calculated finial node
   *               else return 2; not finished
   *
   * @return int
   */
  int state() const;


  /**
   * @brief define the minus operation
   *
   * @param b
   * @return equation<T>
   */
  equation<T> & operator -= (const equation<T> & b);

  /**
     * @brief define the output operation
     *
     * @param eq
     * @return equation<T>
     */
  friend std::ostream& operator << (std::ostream &output,
                                    const equation<T> &eq)
  {
    if(eq.e_vec_.empty()){
      output << "# ------ empty expressions with value = " << eq.value() << std::endl;
    }else{
      output << "# ------ SUM coeff * index , val" << std::endl
             << "# ------ ";
      for(equation<T>::eq_const_iterator eqcit = eq.begin(); eqcit != eq.end();){
        const expression<T> & exp = *eqcit;
        output << exp.coefficient << "*X" << exp.index;
        ++eqcit;
        if(eqcit == eq.end())
          output << " = ";
        else
          output << " + ";
      }
      output << eq.value() << std::endl;
    }
    return output;
  }


  /**
   * @brief get the first expression idx
   *
   * @return size_t
   */
  size_t get_prim_idx() const {
    if(e_vec_.empty())
      return -1;
    else
      return e_vec_.front().index;
  }

  int add_expression(const expression<T> & exp){
    e_vec_.push_back(exp);
    return 0;
  }

  T& value() {return value_;}
  const T& value() const {return value_;}

  std::list<expression<T> > e_vec_;
private:
  T value_;
};

template <typename T>
int equation<T>::merge_similar(){
  typedef typename std::list<expression<T> >::iterator leit;
  leit current = e_vec_.begin();
  leit next = current;
  ++next;
  while(next != e_vec_.end()){
    if(next->index > current->index) {
      current = next++;
      continue;
    }else if(next->index == current->index){
      current->coefficient += next->coefficient;
      e_vec_.erase(next++);
      if(fabs(current->coefficient) < 1e-6){
        e_vec_.erase(current++);
        ++next;
      }
    }else{
      std::cerr << "# [error] merge similar function should only be called "
                << "after sorting." << std::endl;
      return __LINE__;
    }
  }
  return 0;
}

template <typename T>
int equation<T>::normalization()
{
  const T coeff = e_vec_.front().coefficient;
  if(fabs(coeff) < 1e-6){
    if(e_vec_.empty()){
      //std::cerr << "# [info] this equation is empty." << std::endl;
      return 0;
    }else
      std::cerr << "# [error] this expression should be removed." << std::endl;
    return __LINE__;
  }

  value() /= coeff;
  for(typename std::list<expression<T> >::iterator lit = e_vec_.begin();
      lit != e_vec_.end(); ++lit){
    expression<T> & ep = *lit;
    ep.coefficient /= coeff;
  }
  return 0;
}

template <typename T>
int equation<T>::state() const
{
  if(e_vec_.empty()){
    if(fabs(value()) < 1e-8)
      return 0; // is cleared
    return -1; // is conflicted
  }else{
    if(e_vec_.size() == 1)
      return 1; // finial variant
    else
      return 2; // not ready
  }
}

template <typename T>
equation<T> & equation<T>::operator -= (const equation<T> & b )
{
  if(&b == this) {
    e_vec_.clear();
    value() = 0;
    return *this;
  }

  for(typename std::list<expression<T> >::const_iterator lecit_b =
      b.e_vec_.begin(); lecit_b != b.e_vec_.end(); ++lecit_b){
    const expression<T> & exp = *lecit_b;
    const size_t &node_idx = exp.index;
    assert(fabs(exp.coefficient) > 1e-6);
    bool is_found = false;

    for(typename std::list<expression<T> >::iterator leit_a = e_vec_.begin();
        leit_a != e_vec_.end(); ++leit_a){
      expression<T> & exp_a = *leit_a;
      if(exp_a.index == node_idx){
        exp_a.coefficient -= exp.coefficient;
        // zeros
        if(fabs(exp_a.coefficient) < 1e-6) {
          e_vec_.erase(leit_a);
          is_found = true;
          break;
        }
      }
    }
    if(!is_found){
      e_vec_.push_back(make_expression(node_idx, -1 * exp.coefficient));
    }
  }

  value() -= b.value();
  sort_equation();
  normalization();

  return *this;
}

template <typename T>
int equation<T>::update(const size_t & node_idx, const  T & node_value)
{
  int changes = 0;
  for(typename std::list<expression<T> >::iterator leit = e_vec_.begin();
      leit != e_vec_.end();){
    expression<T> & exp = *leit;
    if(exp.index == node_idx){
      value() -= exp.coefficient * node_value;
      e_vec_.erase(leit++);
      ++changes;
    }
    ++leit;
  }
  return changes;
}

//! @brief this class only handle Ai+Bi=Ci
template <typename T>
class gauss_eliminator{
  BOOST_MPL_ASSERT_MSG((boost::is_same<T,double>::value ) ||
                       (boost::is_same<T,float>::value ),
                       NON_FLOAT_TYPES_ARE_NOT_ALLOWED, (void));
public:
  /**
 * @brief construct gauss_eliminator class
 *
 * @param nodes input nodes
 * @param node_flag input node_flag which will be tagged as true if the
 *        corresponding node is known
 */
  gauss_eliminator(std::vector<T> & nodes,
                   boost::dynamic_bitset<> & node_flag)
  //std::vector<bool> & node_flag)
    :nodes_(nodes), node_flag_(node_flag){
    idx2equation_.resize(nodes_.size());
  }

  /**
   * @brief add equation to gauss_eliminator, every time an equation is added,
   *        eliminate function is called.
   *
   * @param input equation
   * @return int
   */
  int add_equation(const equation<T> & e);

  /**
   * @brief This function will start to eliminate equations above all added equations
   *
   * @return int return 0 if works fine, or return non-zeros
   */
  int eliminate();


  /**
   * @brief update the equation, it will check all variant, if a variant is
   *        already known, update this equation
   *
   * @param eq input equation
   * @return int return 0 if nothing changes, or retunr 1;
   */
  int update_equation(equation<T> & eq);

  /**
   * @brief This function is used to check the gauss_elimination is valid or not.
   *        Warning!!! This function costs a lot.
   *
   * @param eq input equation
   * @return int
   */
  bool is_valid()const;

private:
  std::vector<T> & nodes_;
  boost::dynamic_bitset<> & node_flag_;
  std::list<equation<T> > es;

  typedef typename std::list<equation<T> >::iterator equation_ptr;
  std::vector<std::list<equation_ptr> > idx2equation_;

  // this map store the smallest expression
  typedef typename std::map<size_t, std::list<equation_ptr> >::iterator prime_eq_ptr;
  std::map<size_t, std::list<equation_ptr> > prime_idx2equation_;
private:

  int gather_variant_index_of_equation(const equation<T> & eq,
                                       std::vector<size_t> & gather)const;
};

template <typename T>
int gauss_eliminator<T>::add_equation(const equation<T> & e){
  es.push_back(e);
  equation<T> & e_back = es.back();
  for(typename equation<T>::eq_iterator eit = e_back.begin(); eit != e_back.end(); ){
    if(node_flag_[eit->index]){
      e_back.value() -= nodes_[eit->index] * eit->coefficient;
      e_back.e_vec_.erase(eit++);
    }else
      ++eit;
  }
  if(e_back.state() == 0){// this equation is cleared
    es.pop_back();
    return 0;
  }else if(e_back.state() == -1){
    std::cerr << "# [error] strange conflict equation: " << std::endl;
    std::cerr << e;
    es.pop_back();
    return __LINE__;
  }
  e_back.standardization();

  for(typename equation<T>::eq_const_iterator it = e_back.begin();
      it != e_back.end(); ++it){
    equation_ptr end_ptr = es.end();
    idx2equation_[it->index].push_back(--end_ptr);
  }

  equation_ptr end_ptr = es.end();
  prime_idx2equation_[e_back.get_prim_idx()].push_back(--end_ptr);
  eliminate();
  return 0;
}

template <typename T>
int gauss_eliminator<T>::update_equation(equation<T> & eq)
{
  for(typename std::list<expression<T> >::iterator it = eq.e_vec_.begin();
      it != eq.e_vec_.end(); ){
    const expression<T> & exp = *it;
    if(node_flag_[exp.index]){
      eq.value() -= exp.coefficient * nodes_[exp.index];
      eq.e_vec_.erase(it++);
    }else
      ++it;
  }
  eq.standardization();
  return 0;
}

template <typename T>
bool gauss_eliminator<T>::is_valid()const
{
  // check idx2equations
  for(size_t pi = 0; pi < idx2equation_.size(); ++pi){
    const std::list<equation_ptr> & node_link_eq = idx2equation_[pi];
    if(node_link_eq.empty()) continue;
    {
      for(typename std::list<equation_ptr>::const_iterator cit =
          node_link_eq.begin(); cit != node_link_eq.end(); ++cit){
        const equation<T> & eq = *(*cit);
        //size_t ei = 0;
        typename std::list<expression<T> >::const_iterator lecit =
            eq.e_vec_.begin();
        for(; lecit != eq.e_vec_.end(); ++lecit){
          if(lecit->index == pi)
            break;
        }
        if(lecit == eq.e_vec_.end()) {
          std::cerr << "# [error] can not find node " << pi
                    <<  " in its linking equations." << std::endl;
          return false;
        }
      }
    }
  }

  typedef typename std::map<size_t, std::list<equation_ptr> >::const_iterator mslcit;
  for(mslcit it = prime_idx2equation_.begin(); it != prime_idx2equation_.end();
      ++it){
    const std::list<equation_ptr> & node_link_eq = it->second;
    const std::list<equation_ptr> & node_link_in_vec = idx2equation_[it->first];
    const size_t &node_idx = it->first;
    if(node_link_eq.empty()) continue;
    {
      for(typename std::list<equation_ptr>::const_iterator lcit =
          node_link_eq.begin(); lcit != node_link_eq.end(); ++lcit){
        if(find(node_link_in_vec.begin(), node_link_in_vec.end(), *lcit)
           == node_link_in_vec.end()){
          std::cerr << "# [error] can not find node " << node_idx
                    <<  " in its linking equations." << std::endl;
          return false;
        }
      }
    }
  }
  return true;
}


template <typename T>
int gauss_eliminator<T>::gather_variant_index_of_equation(
    const equation<T> & eq, std::vector<size_t> & gather) const
{
  gather.clear();
  for(typename std::list<expression<T> >::const_iterator it = eq.e_vec_.begin();
      it !=  eq.e_vec_.end(); ++it){
    const expression<T> & exp = *it;
    gather.push_back(exp.index);
  }
  return 0;
}

template <typename T>
int gauss_eliminator<T>::eliminate()
{
  while(1){
    bool is_modified = false;

    for(prime_eq_ptr ptr = prime_idx2equation_.begin();
        ptr != prime_idx2equation_.end();)
    {
      std::list<equation_ptr> & dle = ptr->second;
      if(dle.empty()) {
        prime_idx2equation_.erase(ptr++);
        is_modified = true;
        continue;
      }else if(dle.size() == 1){ // contain only one equation
        const int state_ = dle.front()->state();
        if(state_ == 0) { // cleared
          es.erase(dle.front());
          prime_idx2equation_.erase(ptr);
          is_modified = true;
          break;
        }else if(state_ == -1){ // conflict equation
          std::cerr << "# [error] conflict equation " << std::endl;
          return __LINE__;
        }else if(state_ == 1){ // finial variant
          const equation<T> & eq = *dle.front();
          const T &value_ = eq.value();
          // prime index's coefficient should be 1
          assert(fabs(eq.e_vec_.front().coefficient - 1) < 1e-6);
          const size_t index = eq.get_prim_idx();
          if(node_flag_[index]){
            if(fabs(nodes_[index] - value_) > 1e-6){
              std::cerr << "# [error] conficts happen, node " << index
                        << " has different value " << nodes_[index] << ","
                        << value_ << std::endl;
              return __LINE__;
            }
          }else{
            // update corresponding equations with the new node value
            nodes_[index] = value_;
            node_flag_[index] = true;
            std::list<equation_ptr> & node_linked_eq = idx2equation_[index];
            for(typename std::list<equation_ptr>::iterator leqit =
                node_linked_eq.begin(); leqit != node_linked_eq.end();){
              equation<T> & eq = *(*leqit);
              node_linked_eq.erase(leqit++);
              eq.standardization();
            }
            //node_linked_eq.clear();
            prime_idx2equation_.erase(ptr);
            is_modified = true;
          }
        }
        ++ptr;
      }else{
        assert(dle.size() > 1);
        // this prime_index point to several equations,
        // which sould be eliminated
        typename std::list<equation_ptr>::iterator begin = dle.begin();
        typename std::list<equation_ptr>::iterator first = begin++;
        // to keep each prime index linked only one equation
        for(typename std::list<equation_ptr>::iterator next = begin;
            next != dle.end();){

          // gather all index which belong to *next equation,
          // and check if these will change after minus operation
          std::vector<size_t> index_collec_prev;
          gather_variant_index_of_equation(*(*next), index_collec_prev);

          // to eliminate the prime index, each equation minus the first one
          *(*next) -= *(*first);
          (*next)->standardization();

          std::vector<size_t> index_collec_after;
          gather_variant_index_of_equation(*(*next), index_collec_after);

          // if an variant which belongs to next, and vanish after minus operation,
          // should remove the equation linking from idx2equations
          for(size_t idx = 0; idx < index_collec_prev.size(); ++idx){
            if(find(index_collec_after.begin(),
                    index_collec_after.end(),
                    index_collec_prev[idx]) == index_collec_after.end())
            {
              std::list<equation_ptr> & node_linked_eq =
                  idx2equation_[index_collec_prev[idx]];//.erase((*next));
              typename std::list<equation_ptr>::iterator leqit =
                  find(node_linked_eq.begin(),
                       node_linked_eq.end(),(*next));
              if(leqit == node_linked_eq.end()){
                std::cerr << "# [error] can not find this equation of node "
                          << index_collec_prev[idx] << std::endl;
                return __LINE__;
              }else
                node_linked_eq.erase(leqit);
            }
          }

          // if an variant which belong to first, and exist in next after minus
          // operation, should add them into the idx2equations
          for(size_t idx = 0; idx < index_collec_after.size(); ++idx){
            if(find(index_collec_prev.begin(), index_collec_prev.end(),
                    index_collec_after[idx]) == index_collec_prev.end()){
              idx2equation_[index_collec_after[idx]].push_back(*next);
            }
          }

          if((*next)->state() == 0) {// empty equation
            dle.erase(next++);
          }else{
            const size_t prim_index = (*next)->get_prim_idx();
            if(prim_index > (*first)->get_prim_idx()){
              prime_idx2equation_[prim_index].push_back(*next);
              dle.erase(next++);
            }else{
              assert(prim_index == (*first)->get_prim_idx());
              ++next;
            }
          }
        }
        ++ptr;
      } // end else
    }
    if(!is_modified) break;
  }
  return 0;
}
}
}
#endif // GAUSS_ELIMINATION_H
