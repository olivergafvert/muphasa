//
//  landscapes.h
//  mph
//
//  Created by Oliver Gafvert on 04/10/2022.
//  Copyright © 2022 Oliver. All rights reserved.
//

#ifndef landscapes_h
#define landscapes_h

typedef std::vector<std::pair<grade_t, std::vector<size_t>>> MultigradedBasis;
typedef std::vector<std::pair<grade_t, size_t>> Landscape;

class CompressedLandscape : public std::vector<std::pair<grade_t, std::vector<grade_t>>>{
private:
public:
   
    CompressedLandscape() : std::vector<std::pair<grade_t, std::vector<grade_t>>>() { }

    void optimize();
    
    void sort_colex();
    
    size_t operator()(grade_t& grade, size_t k);
};

inline void CompressedLandscape::optimize() {
    for (size_t i=0; i<this->size(); i++) {
        std::vector<grade_t> min_gens;
        std::vector<grade_t>& orig = this->at(i).second;
        for (size_t j=0; j<orig.size(); j++) {
            bool should_add = true;
            for (size_t k=j+1; k<orig.size(); k++) {
                if ( orig[j].leq_poset(orig[k]) ){
                    should_add = false;
                    break;
                }
            }
            if (should_add) {
                min_gens.push_back(orig[j]);
            }
        }
        sort(min_gens.begin(), min_gens.end(), [ ](  grade_t& lhs,  grade_t& rhs )
             {
                 return lhs.lt_colex(rhs) ;
             });
        this->at(i).second = min_gens;
    }
}

inline void CompressedLandscape::sort_colex() {
    for (size_t i=0; i<this->size(); i++) {
        sort(this->at(i).second.begin(), this->at(i).second.end(), [ ](  grade_t& lhs,  grade_t& rhs )
             {
                 return lhs.lt_colex(rhs) ;
             });
    }
}

inline size_t CompressedLandscape::operator()(grade_t& grade, size_t k)
{
    if (k<=0) {
        return 0;
    }
    size_t r_index = grade.size()-1;
    size_t v_r = grade[r_index];
    std::vector<size_t> diffs;
    for(auto& entry : *this) {
        if (entry.first.leq_poset(grade)) {
            size_t low_index = entry.first[r_index];
            grade_t h_grade(grade);
            size_t high_index = low_index;
            if ( entry.second.size() > 0 ){
                bool should_add = false;
                for(auto& _grade : entry.second) {
                    if ( _grade[r_index] >= grade[r_index]) {
                        h_grade[r_index] = _grade[r_index];
                        grade_t g_join = grade.join(_grade);
                        if( g_join.leq_poset(h_grade) ) {
                            high_index = g_join[r_index]-1 >= v_r ? g_join[r_index]-1 : v_r;
                            should_add = true;
                            break;
                        }
                    } else if ( _grade.leq_poset(grade) ) {
                        high_index = v_r;
                        should_add = true;
                        break;
                    }
                }
                if ( should_add && high_index >= v_r ) {
                    diffs.push_back( v_r-low_index > high_index-v_r ? high_index-v_r : v_r-low_index );
                } else {
                    diffs.push_back(v_r-low_index);
                }
            } else {
                diffs.push_back(v_r-low_index);
            }
            
        }
    }
    sort(diffs.begin(), diffs.end());
    if ( diffs.size() > 0 ) {
        std::cout << std::endl;
    }
    return k < diffs.size() ? diffs[diffs.size()-1-k] : 0;
}



#endif /* landscapes_h */
