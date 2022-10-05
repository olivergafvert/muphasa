#ifndef MPH_IO_INCLUDED
#define MPH_IO_INCLUDED

#include <iostream>
#include <vector>
#include "utils.h"
#include "signatureColumn.h"


//TODO: fix the virtual/override of the metrics and filters.
class Metric{
    /** Abstract class representing a metric. **/
public:
    virtual input_t eval(std::vector<input_t> a, std::vector<input_t> b) {return 0;}
};

// TODO: implement more types of metrcs.
class SquaredEuclideanMetric : public Metric{
public:
    virtual input_t eval(std::vector<input_t> a, std::vector<input_t> b) override{
        assert(a.size()==b.size());
        input_t ret = 0;
        for(size_t i=0; i<a.size(); i++){
            ret += (a[i]-b[i])*(a[i]-b[i]);
        }
        return ret;
    }
};


class Filter{
    /** Abstract class representing a filter function.
     A filter function is a function f: R^N -> R taking a point to its filtration value. **/
public:
    virtual input_t eval(std::vector<input_t> a) {return 0;}
};

// TODO: Implement more types of filters.
class XFilter : public Filter{
public:
    int axis;
    input_t coeff;
    XFilter(int _axis, input_t _coeff = 1) : Filter(){
        axis = _axis;
        coeff = _coeff;
    }
    virtual input_t eval(std::vector<input_t> a) override{
        return a[this->axis]*coeff;
    }
};

class DensityFilter : public Filter{
public:
    std::vector<std::vector<input_t>> points;
    hash_map<std::vector<input_t>*, input_t> cache;
    DensityFilter(std::vector<std::vector<input_t>>& _points) : Filter(){
        points = _points;
    }
    virtual input_t eval(std::vector<input_t> a) override{
        //if(cache.find(&a) != cache.end()){
        //    return cache[&a];
        //}
        std::vector<input_t> neighbourlist;
        Metric* m = new SquaredEuclideanMetric();
        for(auto& p:points){
            neighbourlist.push_back(m->eval(std::vector<input_t>(&a[0],&a[p.size()]), p));
        }
        std::sort(neighbourlist.begin(), neighbourlist.end());
        input_t density = 0;
        for(size_t i=0; i<neighbourlist.size(); i++){
            if(neighbourlist[i] > 0){
                input_t d = (i+1.0)/std::pow(neighbourlist[i], points[0].size());
                if(d>density){
                    density = d;
                }
            }
        }
        //cache[&a] = density;
        return density;
    }
};


template <typename T> class Simplex{
public:
    std::vector<size_t> vertices;
    std::vector<T> grade;
    size_t signature = -1;
    
    Simplex(std::vector<size_t> _vertices, std::vector<T> _grade){
        vertices = _vertices;
        grade = _grade;
    }
    
    Simplex(std::vector<size_t> _vertices, std::vector<T> _grade, size_t _signature){
        vertices = _vertices;
        grade = _grade;
        signature = _signature;
    }
};

class VertexContainer : public std::vector<size_t>{
public:
    VertexContainer() : std::vector<size_t>(){}
    VertexContainer(std::vector<size_t>& v) : std::vector<size_t>(v){}
    
    bool operator==(const VertexContainer& other) const {
        for(size_t index=0; index<this->size(); index++){
            if(this->at(index) != other[index]){
                return false;
            }
        }
        return true;
    }
};

struct VectorHasher
{
    /*size_t operator()(std::vector<size_t> const& vec) const {
        size_t seed = vec.size();
        for(auto x : vec) {
            x = ((x >> 16) ^ x) * 0x45d9f3b;
            x = ((x >> 16) ^ x) * 0x45d9f3b;
            x = (x >> 16) ^ x;
            seed ^= x + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }*/
    std::size_t operator()(const VertexContainer& k) const
    {
        return boost::hash_range(k.begin(), k.end());
    }
};

void build_VR_subcomplex(std::vector<std::vector<input_t>>& points, std::vector<Metric*>& metrics, std::vector<Filter*>& filters, std::vector<size_t>& vertices, std::vector<input_t>& prev_values, std::vector<input_t>& max_metric_values, std::vector<Simplex<input_t>>& low_simplices, std::vector<Simplex<input_t>>& mid_simplices, std::vector<Simplex<input_t>>& high_simplices, int hom_dim){
    /* Helper function to recursively compute the Vietoris-Rips complex.
     
     points {std::vector<std::vector<input_t>>} -- a list of points in space.
     metrics {std::vector<Metric>} -- a list of metrics on the points that should be used to contruct the complex.
     filters {std::vector<Filter>} -- a list of filter functions on the points that should be used to contruct the complex.
     vertices {std::vector<size_t>} -- a list of the current vertices in the simplex being contructed.
     prev_values {std::vector<input_t>} -- the grade of the current simplex being contructed.
     max_metric_values {std::vector<input_t>} -- the maximum allowed distances for each metric respectively.
     low_simplices {std::vector<Simplex>} -- all simplices of dimension hom_dim-1.
     mid_simplices {std::vector<Simplex>} -- all simplices of dimension hom_dim.
     high_simplices {std::vector<Simplex>} -- all simplices of dimension hom_dim+1.
     hom_dim {int} -- the dimension of the homology to be computed.
     
     */
    if (vertices.size() == hom_dim) //simplex of dimension hom_dim - 1
    {
        low_simplices.push_back(Simplex<input_t>(vertices, prev_values));
    } else if (vertices.size() == hom_dim + 1) //simplex of dimension hom_dim
    {
        mid_simplices.push_back(Simplex<input_t>(vertices, prev_values));
    } else if (vertices.size() == hom_dim + 2) //simplex of dimension hom_dim + 1
    {
        high_simplices.push_back(Simplex<input_t>(vertices, prev_values));
        return;
    }
    
    for (size_t j = vertices.back() + 1; j < points.size(); j++) {
        std::vector<input_t> curr_values(prev_values);
        bool should_add = true;
        for (size_t k = 0; k < vertices.size(); k++) {
            size_t p = vertices[k];
            for(size_t metric_index = 0; metric_index<metrics.size(); metric_index++){
                input_t d = (*metrics[metric_index]).eval(points[p], points[j]);
                if(d < max_metric_values[metric_index]){
                    curr_values[metric_index] = d < curr_values[metric_index] ? curr_values[metric_index] : d;
                }else{
                    should_add=false;
                    break;
                }
            }
            if(!should_add){
                break;
            }
        }
        
        if(should_add){
            for(size_t filter_index=0; filter_index<filters.size(); filter_index++){
                curr_values[metrics.size()+filter_index] = (*filters[filter_index]).eval(points[j]) < curr_values[metrics.size()+filter_index] ? curr_values[metrics.size()+filter_index] : (*filters[filter_index]).eval(points[j]);
            }
            vertices.push_back(j);
            build_VR_subcomplex(points, metrics, filters, vertices, curr_values, max_metric_values, low_simplices, mid_simplices, high_simplices, hom_dim);
            vertices.pop_back();
        }
    }
}

void build_VR_subcomplex_spatiotemporal(std::vector<std::vector<std::vector<input_t>>>& trajectory_dms, Metric* metric, std::vector<size_t>& vertices, std::vector<input_t>& prev_values, input_t& max_metric_value, std::vector<Simplex<input_t>>& low_simplices, std::vector<Simplex<input_t>>& mid_simplices, std::vector<Simplex<input_t>>& high_simplices, int hom_dim){
    /* Helper function to recursively compute the Vietoris-Rips complex for a dynamic metric space in the sense of Memoli/Kim (2020).
     
     points {std::vector<std::vector<input_t>>} -- a list of points in space.
     metrics {std::vector<Metric>} -- a list of metrics on the points that should be used to contruct the complex.
     filters {std::vector<Filter>} -- a list of filter functions on the points that should be used to contruct the complex.
     vertices {std::vector<size_t>} -- a list of the current vertices in the simplex being contructed.
     prev_values {std::vector<input_t>} -- the grade of the current simplex being contructed.
     max_metric_values {std::vector<input_t>} -- the maximum allowed distances for each metric respectively.
     low_simplices {std::vector<Simplex>} -- all simplices of dimension hom_dim-1.
     mid_simplices {std::vector<Simplex>} -- all simplices of dimension hom_dim.
     high_simplices {std::vector<Simplex>} -- all simplices of dimension hom_dim+1.
     hom_dim {int} -- the dimension of the homology to be computed.
     
     */
    if (vertices.size() == hom_dim) //simplex of dimension hom_dim - 1
    {
        low_simplices.push_back(Simplex<input_t>(vertices, prev_values));
    } else if (vertices.size() == hom_dim + 1) //simplex of dimension hom_dim
    {
        mid_simplices.push_back(Simplex<input_t>(vertices, prev_values));
    } else if (vertices.size() == hom_dim + 2) //simplex of dimension hom_dim + 1
    {
        high_simplices.push_back(Simplex<input_t>(vertices, prev_values));
        return;
    }
    
    for (size_t j = vertices.back() + 1; j < trajectory_dms[0].size(); j++) {
        std::vector<input_t> curr_values(prev_values);
        bool should_add = true;
        
        for(size_t time_index = 0; time_index<trajectory_dms.size(); time_index++){
            should_add = true;
            for (size_t k = 0; k < vertices.size(); k++) {
                size_t p = vertices[k];
                input_t d = trajectory_dms[time_index][j][p];
                if(d < max_metric_value){
                    curr_values[0] = d < curr_values[0] ? curr_values[0] : d;
                }else{
                    should_add=false;
                    break;
                }
            }
            if(should_add){
                break;
            }
        }
        
        if(should_add){
            vertices.push_back(j);
            build_VR_subcomplex_spatiotemporal(trajectory_dms, metric, vertices, curr_values, max_metric_value, low_simplices, mid_simplices, high_simplices, hom_dim);
            vertices.pop_back();
        }
    }
}

void build_VR_subcomplex_grades(std::vector<std::vector<input_t>>& points, std::vector<grade_t>& point_grades, std::vector<size_t>& vertices, grade_t& prev_grade, grade_t& max_grade, std::vector<Simplex<index_t>>& low_simplices, std::vector<Simplex<index_t>>& mid_simplices, std::vector<Simplex<index_t>>& high_simplices, int hom_dim){
    /* Helper function to recursively compute the Vietoris-Rips complex.
     
     points {std::vector<std::vector<input_t>>} -- a list of points in space.
     metrics {std::vector<Metric>} -- a list of metrics on the points that should be used to contruct the complex.
     filters {std::vector<Filter>} -- a list of filter functions on the points that should be used to contruct the complex.
     vertices {std::vector<size_t>} -- a list of the current vertices in the simplex being contructed.
     prev_values {std::vector<input_t>} -- the grade of the current simplex being contructed.
     max_metric_values {std::vector<input_t>} -- the maximum allowed distances for each metric respectively.
     low_simplices {std::vector<Simplex>} -- all simplices of dimension hom_dim-1.
     mid_simplices {std::vector<Simplex>} -- all simplices of dimension hom_dim.
     high_simplices {std::vector<Simplex>} -- all simplices of dimension hom_dim+1.
     hom_dim {int} -- the dimension of the homology to be computed.
     
     */
    if (vertices.size() == hom_dim) //simplex of dimension hom_dim - 1
    {
        low_simplices.push_back(Simplex<index_t>(vertices, prev_grade));
    } else if (vertices.size() == hom_dim + 1) //simplex of dimension hom_dim
    {
        mid_simplices.push_back(Simplex<index_t>(vertices, prev_grade));
    } else if (vertices.size() == hom_dim + 2) //simplex of dimension hom_dim + 1
    {
        high_simplices.push_back(Simplex<index_t>(vertices, prev_grade));
        return;
    }
    
    for (size_t j = vertices.back() + 1; j < points.size(); j++) {
        grade_t grade =  prev_grade.join(point_grades[j]);
        bool should_add = grade.leq_poset(max_grade);
        
        if(should_add){
            vertices.push_back(j);
            build_VR_subcomplex_grades(points, point_grades, vertices, grade, max_grade, low_simplices, mid_simplices, high_simplices, hom_dim);
            vertices.pop_back();
        }
    }
}

void build_VR_subcomplex_dm(std::vector<std::vector<std::vector<input_t>>>& distance_matrices, std::vector<std::vector<input_t>>& filters, std::vector<size_t>& vertices, std::vector<input_t>& prev_values, std::vector<input_t>& max_metric_values, std::vector<Simplex<input_t>>& low_simplices, std::vector<Simplex<input_t>>& mid_simplices, std::vector<Simplex<input_t>>& high_simplices, int hom_dim){
    /* Helper function to recursively compute the Vietoris-Rips complex.
     
     distance_matrices {std::vector<std::vector<std::vector<input_t>>>} -- a list of distance matrices.
     filters {std::vector<std::vector<input_t>>} -- a list of evaluations of filter functions on the points that will be used to contruct the complex.
     vertices {std::vector<size_t>} -- a list of the current vertices in the simplex being contructed.
     prev_values {std::vector<input_t>} -- the grade of the current simplex being contructed.
     max_metric_values {std::vector<input_t>} -- the maximum allowed distances for each metric respectively.
     low_simplices {std::vector<Simplex>} -- all simplices of dimension hom_dim-1.
     mid_simplices {std::vector<Simplex>} -- all simplices of dimension hom_dim.
     high_simplices {std::vector<Simplex>} -- all simplices of dimension hom_dim+1.
     hom_dim {int} -- the dimension of the homology to be computed.
     
     */
    if (vertices.size() == hom_dim) //simplex of dimension hom_dim - 1
    {
        low_simplices.push_back(Simplex<input_t>(vertices, prev_values));
    } else if (vertices.size() == hom_dim + 1) //simplex of dimension hom_dim
    {
        mid_simplices.push_back(Simplex<input_t>(vertices, prev_values));
    } else if (vertices.size() == hom_dim + 2) //simplex of dimension hom_dim + 1
    {
        high_simplices.push_back(Simplex<input_t>(vertices, prev_values));
        return;
    }
    
    size_t n_points = distance_matrices[0].size();
    
    for (size_t j = vertices.back() + 1; j < n_points; j++) {
        std::vector<input_t> curr_values(prev_values);
        bool should_add = true;
        for (size_t k = 0; k < vertices.size(); k++) {
            size_t p = vertices[k];
            for(size_t metric_index = 0; metric_index<distance_matrices.size(); metric_index++){
                input_t d = distance_matrices[metric_index][p][j];
                if(d < max_metric_values[metric_index]){
                    curr_values[metric_index] = d < curr_values[metric_index] ? curr_values[metric_index] : d;
                }else{
                    should_add=false;
                    break;
                }
            }
            if(!should_add){
                break;
            }
        }
        
        if(should_add){
            for(size_t filter_index=0; filter_index<filters.size(); filter_index++){
                curr_values[distance_matrices.size()+filter_index] = filters[filter_index][j] < curr_values[distance_matrices.size()+filter_index] ? curr_values[distance_matrices.size()+filter_index] : filters[filter_index][j];
            }
            vertices.push_back(j);
            build_VR_subcomplex_dm(distance_matrices, filters, vertices, curr_values, max_metric_values, low_simplices, mid_simplices, high_simplices, hom_dim);
            vertices.pop_back();
        }
    }
}

template <typename T> grade_t transform_grade(std::vector<T>& input_grade, std::vector<hash_map<T, index_t>>& grade_map){
    /* Transforms the of the simplex of type 'input_t' to a 'grade_t' type.
     
     Arguments:
        std::vector<T>& input_grade - the input grade to be transformed.
        std::vector<hash_map<T, index_t>>& grade_map - a map mapping 'input_t' type values to the corresponding value in the grade_t type.
     
     Returns:
        grade_t - a transformed grade.
     */
    grade_t grade;
    for(size_t i=0; i<input_grade.size(); i++){
        grade.push_back(grade_map[i][input_grade[i]]);
    }
    return grade;
}

template <typename T> void add_boundary_columns(std::vector<Simplex<T>>& high_simplices, std::vector<Simplex<T>>& low_simplices, Matrix& matrix, std::vector<hash_map<T, index_t>>& grade_map){
    /* Helper function to construct and add the columns of the boundary matrix from lists of simplices. */
    
    std::unordered_map<VertexContainer, size_t, VectorHasher> map;
    for(size_t index=0; index<low_simplices.size(); index++){
        map[VertexContainer(low_simplices[index].vertices)] = index;
    }
    
    for(size_t index=0; index < high_simplices.size(); index++){
        grade_t grade(transform_grade<T>(high_simplices[index].grade, grade_map));
        SignatureColumn column(grade, index);
        sort(high_simplices[index].vertices.begin(), high_simplices[index].vertices.end());
        int coeff = 1;
        for(size_t simplex_index = 0; simplex_index < high_simplices[index].vertices.size() ; simplex_index++){
            VertexContainer vertices;
            std::vector<size_t> high_vertices;
            for(size_t k=0; k < high_simplices[index].vertices.size(); k++){
                if(k!=simplex_index){
                    vertices.push_back(high_simplices[index].vertices[k]);
                    high_vertices.push_back(high_simplices[index].vertices[k]);
                }
            }
            std::vector<size_t> low_vertices = low_simplices[map[vertices]].vertices;
            sort(high_vertices.begin(), high_vertices.end());
            sort(low_vertices.begin(), low_vertices.end());
            if( low_vertices != high_vertices){
                std::cout << "Finished computing boundary matrices." << std::endl;
            }
            column.push(column_entry_t(coeff, map[vertices]));
            //coeff *= -1; //TODO: handle more types of coefficients than F_2.
        }
        matrix.push_back(column);
    }
}

template <typename T> void add_boundary_columns_multicritical(std::vector<Simplex<T>>& high_simplices, std::vector<Simplex<T>>& low_simplices, Matrix& matrix, std::vector<hash_map<T, index_t>>& grade_map){
    /* Helper function to construct and add the columns of the boundary matrix from lists of simplices.
     
     Assumption: low_simplices is ordered lexicographically w.r.t each signature index
     
     */
    
    std::unordered_map<VertexContainer, size_t, VectorHasher> map;
    for(size_t index=0; index<low_simplices.size(); index++){ // Observe that since low_simplices is ordered lex w.r.t each signature index, the representative with lowest lex order will be present in the map.
        map[VertexContainer(low_simplices[index].vertices)] = index;
    }
    
    // Map d
    for(size_t index=0; index < high_simplices.size(); index++){
        grade_t grade(transform_grade<T>(high_simplices[index].grade, grade_map));
        SignatureColumn column(grade, index);
        int coeff = 1;
        for(size_t simplex_index = 0; simplex_index < high_simplices[index].vertices.size() ; simplex_index++){
            VertexContainer vertices;
            std::vector<size_t> high_vertices;
            for(size_t k=0; k < high_simplices[index].vertices.size(); k++){
                if(k!=simplex_index){
                    vertices.push_back(high_simplices[index].vertices[k]);
                    high_vertices.push_back(high_simplices[index].vertices[k]);
                }
            }
            std::vector<size_t> low_vertices = low_simplices[map[vertices]].vertices;
            sort(high_vertices.begin(), high_vertices.end());
            sort(low_vertices.begin(), low_vertices.end());
            if( low_vertices != high_vertices){
                std::cout << "Finished computing boundary matrices." << std::endl;
            }
            size_t diff = 0;
            while ( true ) {
                bool should_add = true;
                for (size_t k=0; k<low_simplices[map[vertices]-diff].grade.size(); k++) {
                    if ( low_simplices[map[vertices]-diff].grade[k]> high_simplices[index].grade[k]) {
                        should_add = false;
                        break;
                    }
                }
                if (should_add) {
                    break;
                } else {
                    diff++;
                }
            }
            size_t m = map[vertices];
            low_vertices = low_simplices[map[vertices]-diff].vertices;
            if( low_vertices != high_vertices){
                throw " Inconsistent basis ";
            }
            column.push(column_entry_t(coeff, map[vertices]-diff));
            //coeff *= -1; //TODO: handle more types of coefficients than F_2.
        }
        matrix.push_back(column);
    }
    
    // Map pi
    std::vector<std::vector<Simplex<T>>> group_by_signature;
    for (size_t index=0; index < low_simplices.size(); index++) {
        size_t signature = low_simplices[index].signature;
        while (group_by_signature.size() <= signature) {
            group_by_signature.push_back(std::vector<Simplex<T>>());
        }
        bool should_add = true;
        for (size_t i=0; i < group_by_signature[signature].size(); i++) {
            bool is_lower = true;
            for (size_t j=0; j < group_by_signature[signature][i].grade.size(); j++) {
                if ( group_by_signature[signature][i].grade[j] > low_simplices[index].grade[j] ) {
                    is_lower = false;
                    break;
                }
            }
            if (is_lower) {
                should_add = false;
                break;
            }
        }
        if (should_add) {
            low_simplices[index].signature = index;
            group_by_signature[signature].push_back(low_simplices[index]);
        }
    }
    
    for (size_t signature=0; signature < group_by_signature.size(); signature++) {
        if (group_by_signature[signature].size() > 1) {
            for (size_t i=0; i < group_by_signature[signature].size(); i++) {
                for (size_t j=0; j < i; j++) {
                    std::vector<T> max_grade;
                    for (size_t t=0; t < group_by_signature[signature][i].grade.size(); t++) {
                        max_grade.push_back( group_by_signature[signature][i].grade[t] < group_by_signature[signature][j].grade[t] ? group_by_signature[signature][j].grade[t] : group_by_signature[signature][i].grade[t] );
                    }
                    grade_t grade(transform_grade<T>(max_grade, grade_map));
                    SignatureColumn column(grade, matrix.size());
                    column.push(column_entry_t(1, group_by_signature[signature][i].signature));
                    column.push(column_entry_t(1, group_by_signature[signature][j].signature));
                    matrix.push_back(column);
                }
            }
        }
    }
    
}


template <typename T> std::vector<hash_map<T, index_t>> compute_grade_map(std::vector<Simplex<T>>& high_simplices, std::vector<Simplex<T>>& mid_simplices, std::vector<Simplex<T>>& low_simplices){
    /* Computes the map needed to transform a 'input_t' type grading of the simplices to the 'grade_t' type used for computations. */

    std::vector<std::set<T>> index_values;
    size_t grade_size = 0;
    if(low_simplices.size()>0){
        grade_size = low_simplices[0].grade.size();
    }
    if(mid_simplices.size() > 0 && mid_simplices[0].grade.size()>grade_size){
        grade_size = mid_simplices[0].grade.size();
    }
    for(size_t i=0; i<grade_size; i++){
        index_values.push_back(std::set<T>());
    }
    
    for (size_t index_column = 0; index_column < high_simplices.size(); index_column++){
        for(size_t index_grade=0; index_grade< high_simplices[index_column].grade.size(); index_grade++){
            index_values[index_grade].insert(high_simplices[index_column].grade[index_grade]);
        }
    }
    for (size_t index_column = 0; index_column < mid_simplices.size(); index_column++){
        for(size_t index_grade=0; index_grade< mid_simplices[index_column].grade.size(); index_grade++){
            index_values[index_grade].insert(mid_simplices[index_column].grade[index_grade]);
        }
    }
    for (size_t index_column = 0; index_column < low_simplices.size(); index_column++){
        for(size_t index_grade=0; index_grade< low_simplices[index_column].grade.size(); index_grade++){
            index_values[index_grade].insert(low_simplices[index_column].grade[index_grade]);
        }
    }
    
    std::vector<std::vector<T>> index_value_lists;
    for(size_t i=0; i<index_values.size(); i++){
        index_value_lists.push_back(std::vector<T>(index_values[i].begin(), index_values[i].end()));
        sort(index_value_lists[i].begin(), index_value_lists[i].end());
    }
    
    std::vector<hash_map<T, index_t>> grade_map;
    for(size_t i=0; i<index_value_lists.size(); i++){
        grade_map.push_back(hash_map<T, index_t>());
        for(size_t j=0; j<index_value_lists[i].size(); j++){
            grade_map[i][index_value_lists[i][j]] = j;
        }
    }
    return grade_map;
}


std::pair<Matrix, Matrix> compute_boundary_matrices(std::vector<std::vector<input_t>>& points, std::vector<Metric*>& metrics, std::vector<Filter*>& filters, std::vector<input_t>& max_metric_values, int hom_dim){
    /* Computes the two boundary matrices of the Vietoris-Rips complex needed to compute the homology of dimension 'hom_dim'.
     
     points {std::vector<std::vector<S>>} -- a list of points in space.
     metrics {std::vector<Metric<S>>} -- a list of metrics on the points that should be used to contruct the complex.
     filters {std::vector<Filter<S>>} -- a list of filter functions on the points that should be used to contruct the complex.
     max_metric_values {std::vector<S>} -- the maximum allowed distances for each metric respectively.
     hom_dim {int} -- the dimension of the homology to be computed.
     
     */
    std::cout << "Starting to compute boundary matrices..." << std::endl;
    Matrix high_matrix;
    Matrix low_matrix;
    
    if(points.size() > 0){
        std::vector<Simplex<input_t>> low_simplices;
        std::vector<Simplex<input_t>> mid_simplices;
        std::vector<Simplex<input_t>> high_simplices;
        
        for (size_t i = 0; i < points.size(); i++) {
            std::vector<size_t> vertices;
            vertices.push_back(i);
            
            std::vector<input_t> curr_values;
            for(size_t metric_index=0; metric_index < metrics.size(); metric_index++){
                curr_values.push_back(0);
            }
            for(size_t filter_index=0; filter_index < filters.size(); filter_index++){
                curr_values.push_back((*filters[filter_index]).eval(points[i]));
            }
            
            std::vector<std::vector<input_t>> points_bak = points;
            
            build_VR_subcomplex(points, metrics, filters, vertices, curr_values, max_metric_values, low_simplices, mid_simplices, high_simplices, hom_dim);
        }
        std::vector<hash_map<input_t, index_t>> grade_map = compute_grade_map<input_t>(high_simplices, mid_simplices, low_simplices);
        add_boundary_columns<input_t>(high_simplices, mid_simplices, high_matrix, grade_map);
        add_boundary_columns<input_t>(mid_simplices, low_simplices, low_matrix, grade_map);
    }
    std::cout << "Finished computing boundary matrices." << std::endl;
    return std::pair<Matrix, Matrix>(high_matrix, low_matrix);
}

std::vector<std::vector<input_t>> compute_distance_matrix(std::vector<std::vector<input_t>>& trajectory, Metric* metric) {
    /* Compressed distance matrix */
    std::vector<std::vector<input_t>> values;
    for (size_t i = 0; i < trajectory.size(); i++) {
        std::vector<input_t> row;
        for (size_t j = 0; j <= i ; j++) {
            input_t val = ((int)(1000*(*metric).eval(trajectory[i], trajectory[j])));
            row.push_back(val);
            std::cout << "Eval: (" << trajectory[i][0] << ", "<< trajectory[i][1] << ", "<< trajectory[i][2] << ")" << "and (" << trajectory[j][0] << ", "<< trajectory[j][1] << ", "<< trajectory[j][2] << ")" << " = " << val << std::endl;
        }
        values.push_back(row);
    }
    return values;
}

std::vector<std::vector<std::vector<input_t>>> get_trajectory_dms(std::vector<std::vector<std::vector<input_t>>>& trajectories, Metric* metric){
    /* Computing the distance matrix at each time-step of the trajectories */
    std::vector<std::vector<std::vector<input_t>>> distance_matrices;
    for (size_t j = 0; j < trajectories[0].size(); j++) {
        std::vector<std::vector<input_t>> trajectory;
        for (size_t i = 0; i < trajectories.size(); i++) {
            trajectory.push_back(trajectories[i][j]);
        }
        distance_matrices.push_back(compute_distance_matrix(trajectory, metric));
    }
    return distance_matrices;
}

std::vector<Simplex<input_t>> calculate_birth_grades2(std::vector<std::vector<std::vector<input_t>>> trajectory_dms, input_t& max_metric_value, std::vector<Simplex<input_t>>& simplices) {
    std::vector<Simplex<input_t>> new_simplices;
    for(size_t index=0; index<simplices.size(); index++){
        input_t prev_value = 0;
        for (size_t i=1; i< simplices[index].vertices.size(); i++) {
            for (size_t j=0; j<i ; j++) {
                prev_value = trajectory_dms[0][i][j] < prev_value ? prev_value : trajectory_dms[0][i][j];
            }
        }
        input_t curr_value = 0;
        for (size_t i=1; i< simplices[index].vertices.size(); i++) {
            for (size_t j=0; j<i ; j++) {
                curr_value = trajectory_dms[1][i][j] < curr_value ? curr_value : trajectory_dms[1][i][j];
            }
        }
        if ( prev_value <= curr_value && prev_value < max_metric_value ) {
            std::vector<input_t> grade;
            grade.push_back(trajectory_dms.size()-1);
            grade.push_back(0);
            grade.push_back(prev_value);
            new_simplices.push_back(Simplex<input_t>(std::vector<size_t>(simplices[index].vertices), std::vector<input_t>(grade), index));
        }
        
        for (size_t time_index=2; time_index < trajectory_dms.size(); time_index++){
            input_t next_value = 0;
            for (size_t i=1; i< simplices[index].vertices.size(); i++) {
                for (size_t j=0; j<i ; j++) {
                    next_value = trajectory_dms[time_index][i][j] < next_value ? next_value : trajectory_dms[time_index][i][j];
                }
            }
            
            if ( curr_value <= prev_value && curr_value <= next_value && curr_value < max_metric_value ) {
                std::vector<input_t> grade;
                grade.push_back(trajectory_dms.size()-time_index);
                grade.push_back(time_index-1);
                grade.push_back(curr_value);
                new_simplices.push_back(Simplex<input_t>(std::vector<size_t>(simplices[index].vertices), std::vector<input_t>(grade), index));
            }
            
            prev_value = curr_value;
            curr_value = next_value;
        }
        
        if ( curr_value <= prev_value && curr_value < max_metric_value ) {
            std::vector<input_t> grade;
            grade.push_back(0);
            grade.push_back(trajectory_dms.size()-1);
            grade.push_back(curr_value);
            new_simplices.push_back(Simplex<input_t>(std::vector<size_t>(simplices[index].vertices), std::vector<input_t>(grade), index));
        }
        
    }
    return new_simplices;
}

std::vector<Simplex<input_t>> calculate_birth_grades(std::vector<std::vector<std::vector<input_t>>> trajectory_dms, input_t& max_metric_value, std::vector<Simplex<input_t>>& simplices) {
    std::vector<Simplex<input_t>> new_simplices;
    for(size_t index=0; index<simplices.size(); index++){
        input_t prev_value = 0;
        for (size_t i=1; i< simplices[index].vertices.size(); i++) {
            for (size_t j=0; j<i ; j++) {
                prev_value = trajectory_dms[0][i][j] < prev_value ? prev_value : trajectory_dms[0][i][j];
            }
        }
        input_t curr_value = 0;
        for (size_t i=1; i< simplices[index].vertices.size(); i++) {
            for (size_t j=0; j<i ; j++) {
                curr_value = trajectory_dms[1][i][j] < curr_value ? curr_value : trajectory_dms[1][i][j];
            }
        }
        if ( prev_value < max_metric_value ) {
            std::vector<input_t> grade;
            grade.push_back(trajectory_dms.size()-1);
            grade.push_back(0);
            grade.push_back(prev_value);
            new_simplices.push_back(Simplex<input_t>(std::vector<size_t>(simplices[index].vertices), std::vector<input_t>(grade), index));
        }
        
        for (size_t time_index=2; time_index < trajectory_dms.size(); time_index++){
            input_t next_value = 0;
            for (size_t i=1; i< simplices[index].vertices.size(); i++) {
                for (size_t j=0; j<i ; j++) {
                    next_value = trajectory_dms[time_index][i][j] < next_value ? next_value : trajectory_dms[time_index][i][j];
                }
            }
            
            if ( curr_value < max_metric_value ) {
                std::vector<input_t> grade;
                grade.push_back(trajectory_dms.size()-time_index);
                grade.push_back(time_index-1);
                grade.push_back(curr_value);
                new_simplices.push_back(Simplex<input_t>(std::vector<size_t>(simplices[index].vertices), std::vector<input_t>(grade), index));
            }
            
            prev_value = curr_value;
            curr_value = next_value;
        }
        
        if ( curr_value < max_metric_value ) {
            std::vector<input_t> grade;
            grade.push_back(0);
            grade.push_back(trajectory_dms.size()-1);
            grade.push_back(curr_value);
            new_simplices.push_back(Simplex<input_t>(std::vector<size_t>(simplices[index].vertices), std::vector<input_t>(grade), index));
        }
        
    }
    return new_simplices;
}

std::pair<Matrix, Matrix> compute_boundary_matrices_spatiotemporal(std::vector<std::vector<std::vector<input_t>>>& trajectories, Metric* metric, input_t max_metric_value, int hom_dim){
    /* Computes the two boundary matrices of the Vietoris-Rips complex needed to compute the homology of dimension 'hom_dim'.
     
     points {std::vector<std::vector<S>>} -- a list of points in space.
     metrics {std::vector<Metric<S>>} -- a list of metrics on the points that should be used to contruct the complex.
     filters {std::vector<Filter<S>>} -- a list of filter functions on the points that should be used to contruct the complex.
     max_metric_values {std::vector<S>} -- the maximum allowed distances for each metric respectively.
     hom_dim {int} -- the dimension of the homology to be computed.
     
     */
    std::cout << "Starting to compute boundary matrices..." << std::endl;
    Matrix high_matrix;
    Matrix low_matrix;
    
    if(trajectories.size() > 0){
        std::vector<Simplex<input_t>> low_simplices;
        std::vector<Simplex<input_t>> mid_simplices;
        std::vector<Simplex<input_t>> high_simplices;
        
        std::vector<std::vector<std::vector<input_t>>> trajectory_dms = get_trajectory_dms(trajectories, metric);
        
        for(size_t i=0; i<trajectory_dms.size(); i++){
            std::cout << "Trajectory dm " << i << ":" << std::endl;
            for (size_t j=0; j<trajectory_dms[i].size(); j++){
                std::cout << "Row " << j << ":" << std::endl;
                for (size_t k=0; k<trajectory_dms[i][j].size(); k++){
                    std::cout << trajectory_dms[i][j][k] << ", ";
                }
                std::cout << std::endl;
            }
        }
        
        for (size_t i = 0; i < trajectories.size(); i++) {
            std::vector<size_t> vertices;
            vertices.push_back(i);
            
            std::vector<input_t> curr_values;
            curr_values.push_back(0);
            
            
            build_VR_subcomplex_spatiotemporal(trajectory_dms, metric, vertices, curr_values, max_metric_value, low_simplices, mid_simplices, high_simplices, hom_dim);
        }
        
        std::vector<Simplex<input_t>> high_simplices_extended = calculate_birth_grades(trajectory_dms, max_metric_value, high_simplices);
        
        std::vector<Simplex<input_t>> mid_simplices_extended = calculate_birth_grades(trajectory_dms, max_metric_value, mid_simplices);
        
        std::vector<hash_map<input_t, index_t>> grade_map = compute_grade_map<input_t>(high_simplices_extended, mid_simplices_extended, low_simplices);
        add_boundary_columns_multicritical<input_t>(high_simplices_extended, mid_simplices_extended, high_matrix, grade_map);
        add_boundary_columns<input_t>(mid_simplices_extended, low_simplices, low_matrix, grade_map);
        std::cout << "Finished computing boundary matrices." << std::endl;
    }
    std::cout << "Finished computing boundary matrices." << std::endl;
    return std::pair<Matrix, Matrix>(high_matrix, low_matrix);
}

std::pair<Matrix, Matrix> compute_boundary_matrices_grades(std::vector<std::vector<input_t>>& points, std::vector<grade_t>& grades, grade_t& max_grade, int hom_dim){
    /* Computes the two boundary matrices of the Vietoris-Rips complex needed to compute the homology of dimension 'hom_dim'.
     
     points {std::vector<std::vector<S>>} -- a list of points in space.
     metrics {std::vector<Metric<S>>} -- a list of metrics on the points that should be used to contruct the complex.
     filters {std::vector<Filter<S>>} -- a list of filter functions on the points that should be used to contruct the complex.
     max_metric_values {std::vector<S>} -- the maximum allowed distances for each metric respectively.
     hom_dim {int} -- the dimension of the homology to be computed.
     
     */
    std::cout << "Starting to compute boundary matrices..." << std::endl;
    Matrix high_matrix;
    Matrix low_matrix;
    
    if(points.size() > 0){
        std::vector<Simplex<index_t>> low_simplices;
        std::vector<Simplex<index_t>> mid_simplices;
        std::vector<Simplex<index_t>> high_simplices;
        
        for (size_t i = 0; i < points.size(); i++) {
            std::vector<size_t> vertices;
            vertices.push_back(i);
            
            grade_t grade = grades[i];
            
            std::vector<std::vector<input_t>> points_bak = points;
            
            build_VR_subcomplex_grades(points, grades, vertices, grade, max_grade, low_simplices, mid_simplices, high_simplices, hom_dim);
        }
        std::vector<hash_map<index_t, index_t>> grade_map = compute_grade_map<index_t>(high_simplices, mid_simplices, low_simplices);
        add_boundary_columns<index_t>(high_simplices, mid_simplices, high_matrix, grade_map);
        add_boundary_columns<index_t>(mid_simplices, low_simplices, low_matrix, grade_map);
    }
    std::cout << "Finished computing boundary matrices." << std::endl;
    return std::pair<Matrix, Matrix>(high_matrix, low_matrix);
}

std::pair<Matrix, Matrix> compute_boundary_matrices_dm(std::vector<std::vector<std::vector<input_t>>>& distance_matrices, std::vector<input_t>& max_metric_values, std::vector<std::vector<input_t>>& filters, int hom_dim){
    /* Computes the two boundary matrices of the Vietoris-Rips complex needed to compute the homology of dimension 'hom_dim'.
     
     distance_matrices {std::vector<std::vector<std::vector<input_t>>>} -- a list of distance matrices.
     max_metric_values {std::vector<S>} -- the maximum allowed distances for each metric respectively.
     filters {std::vector<std::vector<input_t>>} -- a list of evaluations of filter functions on the points that will be used to contruct the complex.
     hom_dim {int} -- the dimension of the homology to be computed.
     
     */
    std::cout << "Starting to compute boundary matrices..." << std::endl;
    Matrix high_matrix;
    Matrix low_matrix;
    
    size_t n_points = distance_matrices[0].size();
    
    if(n_points > 0){
        std::vector<Simplex<input_t>> low_simplices;
        std::vector<Simplex<input_t>> mid_simplices;
        std::vector<Simplex<input_t>> high_simplices;
        
        for (size_t i = 0; i < n_points; i++) {
            std::vector<size_t> vertices;
            vertices.push_back(i);
            
            std::vector<input_t> curr_values;
            for(size_t metric_index=0; metric_index < distance_matrices.size(); metric_index++){
                curr_values.push_back(0);
            }
            for(size_t filter_index=0; filter_index < filters.size(); filter_index++){
                curr_values.push_back(filters[filter_index][i]);
            }
            
            build_VR_subcomplex_dm(distance_matrices, filters, vertices, curr_values, max_metric_values, low_simplices, mid_simplices, high_simplices, hom_dim);
        }
        std::vector<hash_map<input_t, index_t>> grade_map = compute_grade_map<input_t>(high_simplices, mid_simplices, low_simplices);
        add_boundary_columns<input_t>(high_simplices, mid_simplices, high_matrix, grade_map);
        add_boundary_columns<input_t>(mid_simplices, low_simplices, low_matrix, grade_map);
    }
    std::cout << "Finished computing boundary matrices." << std::endl;
    return std::pair<Matrix, Matrix>(high_matrix, low_matrix);
}


void verify_kernel(Matrix& kernel, Matrix& map){
    std::cout << "Starting to verify kernel of map...";
    for(size_t i=0; i<kernel.size(); i++){
        SignatureColumn working_column(kernel[i].grade, i, kernel[i]);
        SignatureColumn result(kernel[i].grade, i);
        while(true){
            index_t pivot = working_column.pop_pivot().get_index();
            if(pivot != -1){
                result.plus(map[pivot]);
            }else{
                break;
            }
        }
        index_t pivot = result.get_pivot().get_index();
        if(pivot!=-1){
            std::cerr << "Column is not in kernel of map";
            throw "Kernel verification failed.";
        }
    }
    std::cout << "Finished verifying kernel of map.";
    std::cout << "Success.";
}

std::vector<std::vector<input_t>> read_point_cloud(std::istream& input_stream, int precision) {
    /** Reads a point cloud from input.
     
     The input should be formatted as: p_11 p_12 ... p_1m
                                       p_21 p_22 ... p_2m etc
     
     Arguments:
            std::istream& input_stream -- input stream of the point cloud
            int precision -- the number of bits of precision TODO: refactor the precision.
     
     Returns:
            std::vector<std::vector<input_t>> -- a list of points represented as vectors with entries of type 'input_t'.
     **/
    std::cout << "Starting to read point cloud...";
    std::vector<std::vector<input_t>> points;
    std::string line;
    input_t value;
    while (std::getline(input_stream, line)) {
        std::vector<input_t> point;
        std::istringstream s(line);
        while (s >> value) {
            point.push_back(input_t(value*pow(10, precision)));
            s.ignore();
        }
        if (!point.empty()){
            points.push_back(point);
        }
    }
    std::cout << "Finished reading points: " << points.size();
    return points;
}

template <typename SignatureColumn, typename Column> void read_input_file(std::istream& input_stream, std::vector<CustomSignatureColumn<SignatureColumn, Column>>& high_matrix, std::vector<CustomSignatureColumn<SignatureColumn, Column>>& low_matrix, int precision=1){
    /**
            Read FIRep formatted file (same as rivet).
     */
    std::string line;
    for(int i=0; i<3; i++){
        std::getline(input_stream, line);
    }
    std::getline(input_stream, line);
    std::istringstream s(line);
    std::vector<int> dims;
    for (int i; s >> i;) {
        dims.push_back(i);
        if (s.peek() == ',' || s.peek() == ' '){
            s.ignore();
        }
    }
    std::vector<std::vector<double>> high_grades;
    std::vector<std::vector<double>> index_value_lists;
    std::vector<hash_map<double, index_t>> visited_map;
    //int precision = 10000000;
    for(int d=0; d<dims[0]; d++){
        std::getline(input_stream, line);
        s = std::istringstream(line);
        bool reading_grade = true;
        CustomSignatureColumn<SignatureColumn, Column> column(grade_t(), d, dims[1], dims[0]);
        high_grades.push_back(std::vector<double>());
        column.syzygy.push(column_entry_t(1, d));
        for (double i; s >> i;) {
            if(reading_grade){
                /*high_grades[d].push_back(i);
                while(index_value_lists.size() < high_grades[d].size()){
                    index_value_lists.push_back(std::vector<double>());
                    visited_map.push_back(hash_map<double, index_t>());
                }
                index_value_lists[high_grades[d].size()-1].push_back(i);*/
                column.grade.push_back((int)(precision*i));
            }else{
                column.push(column_entry_t(1, (int)i));
            }
            while (s.peek() == ' '){
                s.ignore();
            }
            if(s.peek() == ';'){
                s.ignore();
                reading_grade = false;
                while (s.peek() == ' '){
                    s.ignore();
                }
            }
        }
        high_matrix.push_back(column);
    }
    std::vector<std::vector<double>> low_grades;
    for(int d=0; d<dims[1]; d++){
        std::getline(input_stream, line);
        s = std::istringstream(line);
        CustomSignatureColumn<SignatureColumn, Column> column(grade_t(), d, dims[2], dims[1]);
        low_grades.push_back(std::vector<double>());
        column.syzygy.push(column_entry_t(1, d));
        bool reading_grade = true;
        for (double i; s >> i;) {
            if(reading_grade){
                /*low_grades[d].push_back(i);
                while(index_value_lists.size() < low_grades[d].size()){
                    index_value_lists.push_back(std::vector<double>());
                }
                index_value_lists[low_grades[d].size()-1].push_back(i);*/
                column.grade.push_back((int)(precision*i));
            }else{
                column.push(column_entry_t(1, (int)i));
            }
            while (s.peek() == ' '){
                s.ignore();
            }
            if(s.peek() == ';'){
                s.ignore();
                reading_grade = false;
                while (s.peek() == ' '){
                    s.ignore();
                }
            }
        }
        low_matrix.push_back(column);
    }
    
    /*for(size_t i=0; i<index_value_lists.size(); i++){
        sort(index_value_lists[i].begin(), index_value_lists[i].end());
    }
    
    std::vector<hash_map<double, index_t>> grade_map;
    for(size_t i=0; i<index_value_lists.size(); i++){
        grade_map.push_back(hash_map<double, index_t>());
        for(size_t j=0; j<index_value_lists[i].size(); j++){
            grade_map[i][index_value_lists[i][j]] = j;
        }
    }
    for(size_t i=0; i<high_grades.size(); i++){
        for(size_t j=0; j<high_grades[i].size(); j++){
            high_matrix[i].grade.push_back(grade_map[i][high_grades[i][j]]);
        }
    }
    for(size_t i=0; i<low_grades.size(); i++){
        for(size_t j=0; j<low_grades[i].size(); j++){
            low_matrix[i].grade.push_back(grade_map[i][low_grades[i][j]]);
        }
    }*/
}

void apply_diagonal_grade_transform(Matrix columns) {
    for (size_t col_index=0; col_index<columns.size(); col_index++) {
        for (size_t i=0; i<columns[col_index].grade.size()-1; i++) {
            index_t updated_value = columns[col_index].grade[i] - columns[col_index].grade[columns[col_index].grade.size()-1];
            columns[col_index].grade[i] = updated_value >= 0 ? updated_value : 0;
        }
    }
}

#endif // MPH_IO_INCLUDED
