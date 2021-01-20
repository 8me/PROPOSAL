
#include <cmath>
#include <vector>

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Logging.h"
#include "PROPOSAL/geometry/Box.h"

using namespace PROPOSAL;


Box::Box(const Vector3D& position, double x, double y, double z)
    : Geometry("Box", position)
    , x_(x)
    , y_(y)
    , z_(z)
{
}

Box::Box(const nlohmann::json& config)
    : Geometry(config)
{
    if(not config.at("length").is_number())
        throw std::invalid_argument("Length is not a number.");
    if(not config.at("width").is_number())
        throw std::invalid_argument("Width is not a number.");
    if(not config.at("height").is_number())
        throw std::invalid_argument("Height is not a number.");

    x_ = config["length"].get<double>();
    y_ = config["width"].get<double>();
    z_ = config["height"].get<double>();

    if(x_ < 0) throw std::logic_error("lenght must be > 0");
    if(y_ < 0) throw std::logic_error("width must be > 0");
    if(z_ < 0) throw std::logic_error("height must be > 0");
}

bool Box::compare(const Geometry& geometry) const
{
    const Box* box = dynamic_cast<const Box*>(&geometry);

    if (!box)
        return false;
    else if (x_ != box->x_)
        return false;
    else if (y_ != box->y_)
        return false;
    else if (z_ != box->z_)
        return false;
    else
        return true;
}

// ------------------------------------------------------------------------- //
void Box::print(std::ostream& os) const
{
    os << "Width_x: " << x_ << "\tWidth_y " << y_ << "\tHeight: " << z_ << '\n';
}

// ------------------------------------------------------------------------- //
std::pair<double, double> Box::DistanceToBorder(const Vector3D& position, const Vector3D& direction) const
{
    // Calculate intersection of particle trajectory and the box
    // Surface of the box is defined by six planes:
    // E1: x1   =   pos_vec[0] + 0.5*x
    // E2: x1   =   pos_vec[0] - 0.5*x
    // E3: x2   =   pos_vec[1] + 0.5*y
    // E4: x2   =   pos_vec[1] - 0.5*y
    // E5: x3   =   pos_vec[2] + 0.5*z
    // E6: x3   =   pos_vec[2] - 0.5*z
    // straight line (particle trajectory) g = vec(x,y,z) + t * dir_vec( cosph
    // *sinth, sinph *sinth , costh)
    // We are only interested in postive values of t
    // ( we want to find the intersection in direction of the particle
    // trajectory)

    auto dir_vec = direction.GetCartesianCoordinates();
    auto pos_vec = position.GetCartesianCoordinates();

    std::pair<double, double> distance;
    double t;
    double intersection_x;
    double intersection_y;
    double intersection_z;

    std::vector<double> dist;

    double x_calc_pos = position_.GetX() + 0.5 * x_;
    double x_calc_neg = position_.GetX() - 0.5 * x_;
    double y_calc_pos = position_.GetY() + 0.5 * y_;
    double y_calc_neg = position_.GetY() - 0.5 * y_;
    double z_calc_pos = position_.GetZ() + 0.5 * z_;
    double z_calc_neg = position_.GetZ() - 0.5 * z_;

    // intersection with E1
    if (dir_vec[0] != 0) // if dir_vec == 0 particle trajectory is parallel to E1
    {
        t = (x_calc_pos - pos_vec[0]) / dir_vec[0];

        // Computer precision controll
        if (t > 0 && t < GEOMETRY_PRECISION)
            t = 0;

        if (t > 0) // Interection is in particle trajectory direction
        {
            // Check if intersection is inside the box borders
            intersection_y = pos_vec[1] + t * dir_vec[1];
            intersection_z = pos_vec[2] + t * dir_vec[2];
            if (intersection_y >= y_calc_neg && intersection_y <= y_calc_pos && intersection_z >= z_calc_neg &&
                intersection_z <= z_calc_pos)
            {
                dist.push_back(t);
            }
        }
    }

    // intersection with E2
    if (dir_vec[0] != 0) // if dir_vec == 0 particle trajectory is parallel to E2
    {
        t = (x_calc_neg - pos_vec[0]) / dir_vec[0];

        // Computer precision controll
        if (t > 0 && t < GEOMETRY_PRECISION)
            t = 0;

        if (t > 0) // Interection is in particle trajectory direction
        {
            // Check if intersection is inside the box borders
            intersection_y = pos_vec[1] + t * dir_vec[1];
            intersection_z = pos_vec[2] + t * dir_vec[2];
            if (intersection_y >= y_calc_neg && intersection_y <= y_calc_pos && intersection_z >= z_calc_neg &&
                intersection_z <= z_calc_pos)
            {
                dist.push_back(t);
            }
        }
    }

    // intersection with E3
    if (dir_vec[1] != 0) // if dir_vec == 0 particle trajectory is parallel to E3
    {
        t = (y_calc_pos - pos_vec[1]) / dir_vec[1];

        // Computer precision controll
        if (t > 0 && t < GEOMETRY_PRECISION)
            t = 0;

        if (t > 0) // Interection is in particle trajectory direction
        {
            // Check if intersection is inside the box borders
            intersection_x = pos_vec[0] + t * dir_vec[0];
            intersection_z = pos_vec[2] + t * dir_vec[2];
            if (intersection_x >= x_calc_neg && intersection_x <= x_calc_pos && intersection_z >= z_calc_neg &&
                intersection_z <= z_calc_pos)
            {
                dist.push_back(t);
            }
        }
    }

    // intersection with E4
    if (dir_vec[1] != 0) // if dir_vec == 0 particle trajectory is parallel to E4
    {
        t = (y_calc_neg - pos_vec[1]) / dir_vec[1];

        // Computer precision controll
        if (t > 0 && t < GEOMETRY_PRECISION)
            t = 0;

        if (t > 0) // Interection is in particle trajectory direction
        {
            // Check if intersection is inside the box borders
            intersection_x = pos_vec[0] + t * dir_vec[0];
            intersection_z = pos_vec[2] + t * dir_vec[2];
            if (intersection_x >= x_calc_neg && intersection_x <= x_calc_pos && intersection_z >= z_calc_neg &&
                intersection_z <= z_calc_pos)
            {
                dist.push_back(t);
            }
        }
    }

    // intersection with E5
    if (dir_vec[2] != 0) // if dir_vec == 0 particle trajectory is parallel to E5
    {
        t = (z_calc_pos - pos_vec[2]) / dir_vec[2];

        // Computer precision controll
        if (t > 0 && t < GEOMETRY_PRECISION)
            t = 0;

        if (t > 0) // Interection is in particle trajectory direction
        {
            // Check if intersection is inside the box borders
            intersection_x = pos_vec[0] + t * dir_vec[0];
            intersection_y = pos_vec[1] + t * dir_vec[1];
            if (intersection_x >= x_calc_neg && intersection_x <= x_calc_pos && intersection_y >= y_calc_neg &&
                intersection_y <= y_calc_pos)
            {
                dist.push_back(t);
            }
        }
    }

    // intersection with E6
    if (dir_vec[2] != 0) // if dir_vec == 0 particle trajectory is parallel to E6
    {
        t = (z_calc_neg - pos_vec[2]) / dir_vec[2];

        // Computer precision controll
        if (t > 0 && t < GEOMETRY_PRECISION)
            t = 0;

        if (t > 0) // Interection is in particle trajectory direction
        {
            // Check if intersection is inside the box borders
            intersection_x = pos_vec[0] + t * dir_vec[0];
            intersection_y = pos_vec[1] + t * dir_vec[1];
            if (intersection_x >= x_calc_neg && intersection_x <= x_calc_pos && intersection_y >= y_calc_neg &&
                intersection_y <= y_calc_pos)
            {
                dist.push_back(t);
            }
        }
    }

    if (dist.size() < 1) // No intersection with the box
    {
        distance.first  = -1;
        distance.second = -1;
    } else if (dist.size() == 1) // Particle is inside the box and we have one
                                 // intersection in direction of the particle
                                 // trajectory
    {
        distance.first  = dist.at(0);
        distance.second = -1;
    } else if (dist.size() == 2) // Particle is outside and the box is infront
                                 // of the particle trajectory ( two
                                 // intersections).
    {
        if (dist.at(0) == dist.at(1)) {
            // Both distances are equal, this means he have hit an edge of the box
            distance.first = dist.at(0);
            distance.second = -1;
        }
        distance.first  = dist.at(0);
        distance.second = dist.at(1);
        if (distance.second < distance.first)
        {
            std::swap(distance.first, distance.second);
        }

    } else if (dist.size() == 3 && dist.at(0) == dist.at(1) && dist.at(1) == dist.at(2))
    {
        // We have three intersections, all of them are equal, this means we
        // hit an corner of the box
        distance.first = dist.at(0);
        distance.second = -1;
    } else
    {
        Logging::Get("proposal.geometry")->error("This point should never be reached... (-1/-1) is returned");

        distance.first  = -1;
        distance.second = -1;
    }

    // Make a computer precision controll!
    // This is necessary cause due to numerical effects it meight be happen
    // that a particle which is located on a gemoetry border is treated as
    // inside
    // or outside
    if (distance.first < GEOMETRY_PRECISION)
        distance.first = -1;
    if (distance.second < GEOMETRY_PRECISION)
        distance.second = -1;
    if (distance.first < 0)
        std::swap(distance.first, distance.second);

    return distance;
}
