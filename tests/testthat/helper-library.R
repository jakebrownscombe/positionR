# An empty library folder under the session's temporary folder, which R
# removes when the session ends.
empty_library <- function() {
  path <- tempfile("emptylib")
  dir.create(path)
  path
}
