package pwr.adamzimnyy;

/**
 * Created by adamz on 03.06.2017.
 */
public class Rating {
String name;
    int trait;
    double rating;
    int priority;
    @Override
    public String toString() {
        return   "trait=" + trait +
                " rating=" + rating +
                " priority=" + priority+"\n";
    }
}
